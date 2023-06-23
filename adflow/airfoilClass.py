import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.interpolate import splrep, splev
import argparse

parser = argparse.ArgumentParser()
# 1 is x, 2 is y, 3 is z, -1 is -x, -2 is -y, -3 is -z.
parser.add_argument("--chordIndex", type=int, default=1)
parser.add_argument("--liftIndex", type=int, default=2)
parser.add_argument("--spanIndex", type=int, default=3)
parser.add_argument("--sliceFile", type=str, default="input_files/slices.dat")
parser.add_argument("--plotly", type=bool, default=True)
args = parser.parse_args()


def readSlices(filename, nmax=-1):
    """
    Reads slice file output by ADflow and return data and connectivity
    """
    # get all lines
    with open(filename, "r") as f:
        lines = f.readlines()

    # get the number of slices
    n_slice = 0
    for line in lines:
        if "Zone T" in line:
            n_slice += 1

    # max number of slices to read
    if nmax == -1:
        nmax = n_slice

    print(f"Found {n_slice} slices in {filename}")

    # first line is the title
    # second line has the variable list
    var_list = lines[1].replace('"', "").split()
    var_list.remove("Variables")
    var_list.remove("=")
    print("Found variables:", var_list)

    # number of variables
    n_var = len(var_list)

    # list to hold the data and conn arrays for each slice
    slice_data = []
    slice_conn = []
    normals = []

    # beginning line index for the current slice
    slice_begin = 2

    # loop over slices
    for islice in range(nmax):
        slice_header = lines[slice_begin].replace('"', "").replace("(", "").replace(")", "").replace(",", "").split()
        # y_loc = float(slice_header[-1])
        try:
            # if slice is arbitrary, gets the normal
            normal = [float(slice_header[-3]), float(slice_header[-2]), float(slice_header[-1])]
        except Exception:
            # if slice is not arbitrary, normal of slice is spanwise axis
            normal = [0, 0, 0]
            normal[args.spanIndex - 1] = 1
        normals.append(normal)
        slice_counts = lines[slice_begin + 1].replace("=", "").split()
        n_node = int(slice_counts[1])
        n_elem = int(slice_counts[3])
        # print(f"reading slice {islice} at y = {y_loc} with {n_node} nodes and {n_elem} elements")

        # indices where the node data and conn starts
        data_beg = slice_begin + 3
        conn_beg = data_beg + n_node

        # load the node data
        tmp = np.genfromtxt(
            filename,
            dtype={
                "names": var_list,
                "formats": ["f4"] * n_var,
            },
            skip_header=data_beg,
            max_rows=n_node,
        )
        slice_data.append(tmp)
        # print(var_list)
        # load the connectivity
        conn1, conn2 = np.genfromtxt(
            filename,
            dtype={
                "names": (
                    "conn1",
                    "conn2",
                ),
                "formats": ("i4", "i4"),
            },
            usecols=(0, 1),
            skip_header=conn_beg,
            max_rows=n_elem,
            unpack=True,
            invalid_raise=False,
        )

        # -1 is to get back to numpy indexing
        conn = np.column_stack([conn1, conn2]) - 1

        # the line segments in newConn is likely unordered, we need to re-order them
        slice_conn.append(conn)
        slice_begin += 1 + n_node + n_elem + 2
    return slice_data, slice_conn, normals


class Airfoil:
    """
    Create the Airfoil object.

    Parameters
    ----------
    name : str
        name of the slice
    data : dict
        dict output by the cut
    conn : array Nx2
        describe connectivity between points
    normal : array size 3
        normal of the cut
    """

    def __init__(self, name, data, conn, normal):
        self.name = name
        self.data = data
        self.conn = conn
        self.normal = normal
        self.coords = np.zeros((len(self.data["CoordinateX"]), 3))
        self.coords[:, 0] = self.data["CoordinateX"]
        self.coords[:, 1] = self.data["CoordinateY"]
        self.coords[:, 2] = self.data["CoordinateZ"]
        self.te = np.zeros(1)
        self.le = np.zeros(1)
        self.sortConn()
        self.detectTEAngle()
        self.detectLe()
        self.computeLe()

    def FEsort(self):
        """
        This function can be used to sort connectivities coming from the CGNS file

        barsConn should be a list of integers [[2,3],[5,6],[3,4],[5,4],...]

        newConn will be set to None if the sorting algorithm fails.
        """
        # We will solve this in an iterative process until the connectivities do not change
        # anymore
        barsConn = self.conn.tolist()
        # Initialize the newest connectivities as our current input.
        newConn = barsConn

        # Say that we still need to search
        keep_searching = True

        # Get the number of bars
        nBars = len(barsConn)

        # Initialize mapping that will link indices of the given FE data to the sorted one.
        # This mapping will tell which curve the bar element was assigned to.
        # For instance, newMap[5] gives the barIDs that belong to curve[5].
        newMap = [[i] for i in range(nBars)]

        while keep_searching:
            # If nothing happens, we will get out of the loop
            keep_searching = False

            # Update "old" connectivities and mapping
            oldConn = newConn[:]
            oldMap = newMap[:]

            # Initialize list of new connectivities and mapping with the first element of the old ones.
            # We will also pop this element from oldConn. Remember that in subsequent iterations, each
            # element of oldConn represents in fact a concatenation of bar elements.
            newConn = [oldConn.pop(0)]
            newMap = [oldMap.pop(0)]

            # We will keep sorting until all elements of oldConn are assigned and popped out
            # NOTE: An "Element" here may be a concatenation of several bar elements

            while len(oldConn) > 0:
                # Pop another element from oldConn
                oldElement = oldConn.pop(0)
                oldElemMap = oldMap.pop(0)

                # We do not know if this element could be linked beforehand
                linked_element = False

                # Now we check if we can fit this element into any new connectivity
                newElemCounter = 0
                for newElement, newElemMap in zip(newConn, newMap):
                    if oldElement[0] == newElement[0]:
                        # We will flip the old element and place it at the beginning of the
                        # new connectivity
                        newConn[newElemCounter] = oldElement[::-1] + newElement[1:]

                        # Also assign ID the beginning of the map
                        newMap[newElemCounter] = oldElemMap[::-1] + newElemMap[:]

                        # We need to keep searching as we are still making changes
                        linked_element = True

                        # We only need to keep searching if we have updates
                        keep_searching = True

                        # Jump out of the for loop
                        break

                    elif oldElement[0] == newElement[-1]:
                        # Just append the old element
                        newConn[newElemCounter] = newElement[:-1] + oldElement

                        # Also assign mapping
                        newMap[newElemCounter] = newElemMap[:] + oldElemMap

                        # We need to keep searching as we are still making changes
                        linked_element = True

                        # We only need to keep searching if we have updates
                        keep_searching = True

                        # Jump out of the for loop
                        break

                    elif oldElement[-1] == newElement[0]:
                        # Place the old element at the beginning
                        newConn[newElemCounter] = oldElement + newElement[1:]

                        # Also assign mapping
                        newMap[newElemCounter] = oldElemMap + newElemMap[:]

                        # We need to keep searching as we are still making changes
                        linked_element = True

                        # We only need to keep searching if we have updates
                        keep_searching = True

                        # Jump out of the for loop
                        break

                    elif oldElement[-1] == newElement[-1]:
                        # Flip and append the old element
                        newConn[newElemCounter] = newElement[:-1] + oldElement[::-1]

                        # Also assign mapping
                        newMap[newElemCounter] = newElemMap[:] + oldElemMap[::-1]

                        # We need to keep searching as we are still making changes
                        linked_element = True

                        # We only need to keep searching if we have updates
                        keep_searching = True

                        # Jump out of the for loop
                        break

                    # Increment element counter
                    newElemCounter = newElemCounter + 1

                if not linked_element:  # We have an element that is not connected to anyone
                    # Define a new curve in newConn
                    newConn.append(oldElement)

                    # Define a new mapping in newMap
                    newMap.append(oldElemMap)

        # Right now, newConn represent each line by an 1D array of point IDs.
        # e.g. [[2,3,4,5],[1,7,8]]
        # We need to convert this back to FE format, represented by 2D arrays.
        # e.g. [[[2,3,4],[3,4,5]],[[1,7],[7,8]]]

        # Initialize FE array for each curve
        newConnFE = [np.zeros((len(curve) - 1, 2), dtype=int) for curve in newConn]

        # Fill FE data for each curve
        for curveID in range(len(newConn)):
            # Get current curve
            curve = newConn[curveID]

            # Get FE list of the current curve
            FEcurve = newConnFE[curveID]

            # Assign FE connectivity
            for pointID in range(1, len(curve)):
                # Get point indices
                prevPoint = curve[pointID - 1]
                currPoint = curve[pointID]

                # Assign bar FE
                FEcurve[pointID - 1, 0] = prevPoint
                FEcurve[pointID - 1, 1] = currPoint

        # We just need to crop the last index of the mapping arrays because it has a repeated
        # number

        # Now we do a final check to remove degenerate bars (bars that begin and end at the same point)
        # Take every disconnect curve in newConnFE:
        for curveID in range(len(newConnFE)):
            # We convert the connectivity array to list so we can 'pop' elements
            curveFE = newConnFE[curveID].tolist()

            # Get the mapping as well
            currMap = newMap[curveID]

            # Define FE counter
            FEcounter = 0

            while FEcounter < len(curveFE):
                # Get current finite element
                FE = curveFE[FEcounter]

                # Check if the start and end points are the same
                if FE[0] == FE[1]:
                    # Remove FE
                    curveFE.pop(FEcounter)

                    # Remove mapping (remove bar element from this curve)
                    currMap.pop(FEcounter)

                else:
                    # Increment counter
                    FEcounter = FEcounter + 1

            # Convert connectivity back to numpy array
            newConnFE[curveID] = np.array(curveFE)

        # Now remove empty curves
        curveID = 0
        while curveID < len(newConnFE):
            # Check if current curve has no element
            if len(newConnFE[curveID].tolist()) == 0:
                # If this is the case, remove the curve and its mapping
                newConnFE.pop(curveID)
                newMap.pop(curveID)

            else:
                # Increment counter
                curveID = curveID + 1

        # Return the sorted array and the mapping
        return newConnFE, newMap

    def sortConn(self):
        newConn, newMap = self.FEsort()

        # the line segments in newConn is likely unordered, we need to re-order them
        x = self.data["CoordinateX"]
        y = self.data["CoordinateY"]
        z = self.data["CoordinateZ"]

        line_beg = np.zeros((len(newConn), 3))
        line_end = np.zeros((len(newConn), 3))

        for ii, line in enumerate(newConn):
            line_beg[ii, 0] = x[line[0, 0]]
            line_beg[ii, 1] = y[line[0, 0]]
            line_beg[ii, 2] = z[line[0, 0]]
            line_end[ii, 0] = x[line[-1, 1]]
            line_end[ii, 1] = y[line[-1, 1]]
            line_end[ii, 2] = z[line[-1, 1]]

        # now we need to match which line beg node matches which line end node
        # we can just do an np search, since this split should not be that many
        # we start with the first line and match its end node with another line's begin node until we are back at the beginning
        keep_searching = True
        cur_line_ind = 0
        conn_order_map = [0]  # start with the first line so we have the first zero
        while keep_searching:
            cur_line_end = line_end[cur_line_ind]

            # find the closest begin node for the next line (should be an identical node)
            delta = np.linalg.norm(cur_line_end - line_beg, axis=1)
            new_line_ind = np.argmin(delta)

            # check if we have a match. if not, either a line needs to be flipped or this is an open curve
            if np.min(delta) > 1e-12:
                # the next line might need to be flipped. check the line_end array if we get any matches here
                arr = np.linalg.norm(cur_line_end - line_end, axis=1)
                # this should get 2 matches
                matched_inds = np.array(np.where(arr < 1e-12)).squeeze()

                # only keep going if we had 2 matches; one is the current line, the other is the next one
                # TODO check this, was a bug?
                if len(matched_inds) == 2:
                    # this is the index of the next line, but it needs to be flipped
                    new_line_ind = matched_inds[matched_inds != cur_line_ind][0]

                    # flip the beginning and end points
                    tmp = line_end[new_line_ind, :].copy()
                    line_end[new_line_ind, :] = line_beg[new_line_ind, :].copy()
                    line_beg[new_line_ind, :] = tmp

                    # also flip the connectivity array
                    conn_to_flip = np.array(newConn[new_line_ind])

                    # this should flip on all axes
                    flipped_conn = np.flip(conn_to_flip)

                    # put it back
                    newConn[new_line_ind] = flipped_conn

                else:
                    raise Exception(
                        "The distance between the end of the current line and beginning of the next line is larger than 1e-12. I also could not find any line that ends at the same location, so the curve might be open?"
                    )

            # get the new line's end
            new_end = line_end[new_line_ind]

            # append this line to the list
            conn_order_map.append(new_line_ind)

            if np.linalg.norm(new_end - line_beg[0]) < 1e-12:
                # the current line's end is the beginning of the original, so we are done!
                keep_searching = False

            # switch indices
            cur_line_ind = new_line_ind

        # remove the duplicate lines and get an order of nodes. the first and last nodes will be duplicated

        # now append the elements back to back
        final_conn = []
        for ii in range(len(conn_order_map)):
            # extend with the current connectivity
            final_conn.extend(newConn[conn_order_map[ii]].tolist())
        self.sortedConn = np.array(final_conn)

    def detectTEAngle(self):
        """
        Detect the upper, lower trailing edge of the point cloud
        Compute the true Trailing edge as the middle point of the upper, lower TE
        """
        c = self.sortedConn
        v1 = self.coords[c[:, 0], :] - self.coords[c[:, 1], :]
        v1n = np.linalg.norm(v1, axis=1)
        v2 = np.roll(v1, 1, axis=0)
        v2n = np.roll(v1n, 1)
        # find the angles between adjacent vectors
        thetas = np.arccos(np.minimum(np.ones(len(v1n)), np.einsum("ij,ij->i", v1, v2) / (v1n[:] * v2n[:])))
        ids = np.argwhere(thetas > 40.0 * np.pi / 180.0)
        if len(ids) == 0:
            print("WARNING, no sharp corner detected")
        elif len(ids) == 1:
            ind0 = ids[0][0]
            node0 = c[ind0, 0]
            self.lower_te_ind = node0
            self.upper_te_ind = node0
        elif len(ids) == 2:
            ind0 = ids[0][0]
            ind1 = ids[1][0]
            node0 = c[ind0, 0]
            node1 = c[ind1, 0]
            if self.coords[node0, args.liftIndex - 1] > self.coords[node1, args.liftIndex - 1]:
                self.upper_te_ind = node0
                self.lower_te_ind = node1
            else:
                self.lower_te_ind = node0
                self.upper_te_ind = node1
        else:
            print("WARNING, there are more than 2 sharp corners on slice")
        self.upperTe = self.coords[self.upper_te_ind, :]
        self.lowerTe = self.coords[self.lower_te_ind, :]
        self.te = 0.5 * (self.upperTe + self.lowerTe)

    def detectLe(self):
        """
        Detect the leading edge of the existing point cloud as the max distance from the true TE
        """
        c = self.sortedConn[:, 0].tolist()
        te_to_pt_vec = self.coords[c, :] - self.te[:]
        # now compute the distances
        distanmces_to_te = np.linalg.norm(te_to_pt_vec, axis=1)
        # get the max-distance location
        max_dist_ind = np.argmax(distanmces_to_te)
        # save le
        self.le_conn_ind = max_dist_ind
        self.le_ind = c[max_dist_ind]
        self.le = self.coords[self.le_ind, :]
        # indices of nodes from LE to upper blunt Trailing edge
        # indices of nodes from LE to lower blunt Trailing edge
        self.upper_skin_nodes = []
        self.lower_skin_nodes = []

    def reorderData(self):
        """
        Sort data in self.data based on the new connectivity
        """
        conn = self.sortedConn

        # get the node list
        node_list = conn[:, 0].tolist()

        orderdata = {
            "x": self.data["CoordinateX"][node_list],
            "y": self.data["CoordinateY"][node_list],
            "z": self.data["CoordinateZ"][node_list],
            # "u": self.data["VelocityX"][node_list],
            # "v": self.data["VelocityY"][node_list],
            # "w": self.data["VelocityZ"][node_list],
            "cp": self.data["CoefPressure"][node_list],
        }
        self.order_data = orderdata
        # find the LE x index

        # roll all arrays to have LE be the first entry
        for k, v in self.order_data.items():
            self.order_data[k] = np.roll(v, -self.le_conn_ind)

        # check if the orientation is right, we want the upper skin first (just an arbitrary convention)
        # TODO Careful this doesn't work if dihedral over 90 degree
        direction = ["x", "y", "z"]
        if self.order_data[direction[args.liftIndex - 1]][1] - self.order_data[direction[args.liftIndex - 1]][0] > 0:
            # we are going up in the LE, this is what we want
            pass
        else:
            for k, v in self.order_data.items():
                # flip and roll by 1
                self.order_data[k] = np.roll(np.flip(v), 1)

    def computeLe(self):
        """
        For postprocessing twist, chord, etc, compute accurate true LE
        """
        self.reorderData()
        # the LE node and the n-neighboring nodes are interpolated regions with a B-spline
        n_neigh = 10  # number of neighboring nodes to include in both directions

        # fit curve to it
        t = np.linspace(0, 1, num=n_neigh * 2 + 1)
        # construct spline for all data
        tcks = {}
        for k, v in self.order_data.items():
            tcks[k] = splrep(t, np.concatenate((self.order_data[k][-n_neigh:], self.order_data[k][: n_neigh + 1])))

        def min_dist_from_te(t):
            vec = np.array(
                [splev(t, tcks["x"]) - self.te[0], splev(t, tcks["y"]) - self.te[1], splev(t, tcks["z"]) - self.te[2]]
            )
            return -np.linalg.norm(vec)

        # can be improved by providing analytic or CS derivatives
        opt_res = minimize(min_dist_from_te, [0], bounds=[[0, 1]], options={"disp": False}, tol=1e-10)

        if not opt_res["success"]:
            print("WARNING, Optimization failed to find the max distance pt")
        # take the result
        t_opt = opt_res["x"]

        # overwrite the current LE coordinates with the new result
        for k, v in self.order_data.items():
            self.order_data[k][0] = splev(t_opt, tcks[k])
        self.accurate_le = np.array([self.order_data["x"][0], self.order_data["y"][0], self.order_data["z"][0]])

    def getTwist(self):
        untwisted = np.array([0, 0, 0])
        untwisted[int(abs(args.chordIndex)) - 1] = np.sign(args.chordIndex) * 1
        normal = self.normal / np.linalg.norm(self.normal)
        twisted = self.te - self.accurate_le
        twisted = twisted / np.linalg.norm(twisted)
        twist = np.arcsin(np.dot(np.cross(twisted, untwisted), normal))
        self.twist = twist * 180 / np.pi

    def getDihedral(self):
        pass

    def getChord(self):
        if len(self.te) == 3:
            self.chord = np.linalg.norm(self.accurate_le - self.te)
            return self.chord
        else:
            print("WARNING, no le and te detected")
            return 1

    def plotAirfoil(self):
        fig = go.Figure()
        fig.add_trace(
            go.Scatter3d(
                x=self.coords[:, 0],
                y=self.coords[:, 1],
                z=self.coords[:, 2],
                mode="markers",
                marker=dict(color="black", size=4),
            )
        )
        fig.add_trace(
            go.Scatter3d(
                x=[self.lowerTe[0]],
                y=[self.lowerTe[1]],
                z=[self.lowerTe[2]],
                mode="markers",
                marker=dict(color="cyan", size=6),
            )
        )
        fig.add_trace(
            go.Scatter3d(
                x=[self.upperTe[0]],
                y=[self.upperTe[1]],
                z=[self.upperTe[2]],
                mode="markers",
                marker=dict(color="orange", size=6),
            )
        )
        fig.add_trace(
            go.Scatter3d(
                x=[self.accurate_le[0]],
                y=[self.accurate_le[1]],
                z=[self.accurate_le[2]],
                mode="markers",
                marker=dict(color="grey", size=8),
            )
        )
        if len(self.te) == 3:
            fig.add_trace(
                go.Scatter3d(
                    x=[self.te[0]], y=[self.te[1]], z=[self.te[2]], mode="markers", marker=dict(color="blue", size=6)
                )
            )
        if len(self.le) == 3:
            fig.add_trace(
                go.Scatter3d(
                    x=[self.le[0]], y=[self.le[1]], z=[self.le[2]], mode="markers", marker=dict(color="red", size=6)
                )
            )
        fig.update_scenes(aspectmode="data")
        fig.show()


class Wing:
    def __init__(self, list_data, list_conn, list_normal):
        self.datas = list_data
        self.conns = list_conn
        self.normals = list_normal
        self.airfoils = []
        for k in range(len(list_data)):
            airfoil = Airfoil(str(k), list_data[k], list_conn[k], list_normal[k])
            self.airfoils.append(airfoil)

    def getTwistDistribution(self):
        span, twist = [], []
        for air in self.airfoils:
            air.getTwist()
            span.append(air.accurate_le[args.spanIndex - 1])
            twist.append(air.twist)
        self.span = span
        self.twist = twist

    def getChordDistribution(self):
        span, chord = [], []
        for air in self.airfoils:
            air.getChord()
            span.append(air.accurate_le[args.spanIndex - 1])
            chord.append(air.chord)
        self.span = span
        self.chord = chord

    def plotTwistDsitribution(self):
        self.getTwistDistribution()
        plt.plot(self.span, self.twist)
        plt.xlabel("Span (m)")
        plt.ylabel("Twist")
        plt.show()

    def plotChordDsitribution(self):
        wing.getChordDistribution()
        plt.plot(self.span, self.chord)
        plt.xlabel("Span (m)")
        plt.ylabel("Chord (m)")
        plt.show()

    def plotWing(self):
        fig = go.Figure()
        for air in self.airfoils:
            fig.add_trace(
                go.Scatter3d(
                    x=air.coords[:, 0],
                    y=air.coords[:, 1],
                    z=air.coords[:, 2],
                    mode="markers",
                    marker=dict(color="black", size=4),
                )
            )
            fig.add_trace(
                go.Scatter3d(
                    x=[air.lowerTe[0]],
                    y=[air.lowerTe[1]],
                    z=[air.lowerTe[2]],
                    mode="markers",
                    marker=dict(color="cyan", size=6),
                )
            )
            fig.add_trace(
                go.Scatter3d(
                    x=[air.upperTe[0]],
                    y=[air.upperTe[1]],
                    z=[air.upperTe[2]],
                    mode="markers",
                    marker=dict(color="orange", size=6),
                )
            )
            fig.add_trace(
                go.Scatter3d(
                    x=[air.te[0]], y=[air.te[1]], z=[air.te[2]], mode="markers", marker=dict(color="blue", size=6)
                )
            )
            fig.add_trace(
                go.Scatter3d(
                    x=[air.le[0]], y=[air.le[1]], z=[air.le[2]], mode="markers", marker=dict(color="red", size=6)
                )
            )
        fig.update_scenes(aspectmode="data")
        fig.show()


slice_data, slice_conn, normals = readSlices(args.sliceFile)
air = Airfoil("1", slice_data[0], slice_conn[0], normals[0])

wing = Wing(slice_data, slice_conn, normals)
wing.plotTwistDsitribution()
wing.plotChordDsitribution()

if args.plotly:
    import plotly.graph_objects as go

    air.plotAirfoil()
    wing.plotWing()
