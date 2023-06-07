import numpy as np
import warnings
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import plotly.graph_objects as go
from scipy.interpolate import splrep, splev


def readSlices(filename, nmax=-1):
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

    # beginning line index for the current slice
    slice_begin = 2

    # loop over slices
    for islice in range(nmax):
        slice_header = lines[slice_begin].replace('"', "").replace("(", "").replace(")", "").split()
        y_loc = float(slice_header[-1])

        slice_counts = lines[slice_begin + 1].replace("=", "").split()
        n_node = int(slice_counts[1])
        n_elem = int(slice_counts[3])
        print(f"reading slice {islice} at y = {y_loc} with {n_node} nodes and {n_elem} elements")

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
        print(var_list)
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
    return slice_data, slice_conn


slice_data, slice_conn = readSlices("input_files/slice.dat")


class Airfoil:
    """
    Create the Airfoil object.

    Parameters
    ----------
    data : dict
        dict output by the cut
    conn : list Nx2
        describe connectivity between points
    """

    def __init__(self, data, conn):
        self.data = data
        self.conn = conn
        self.coords = np.zeros((len(data["CoordinateX"]), 3))
        self.coords[:, 0] = data["CoordinateX"]
        self.coords[:, 1] = data["CoordinateY"]
        self.coords[:, 2] = data["CoordinateZ"]
        self.val = {} 
        self.val["cp"] = data["CoefPressure"]
        self.val["u"] = data['VelocityX']
        self.val["v"] = data['VelocityY']
        self.val["w"] = data['VelocityZ']

    def sortConn(self):
        """
        This function can be used to sort connectivities coming from the CGNS file
        barsConn should be a list of integers [[2,3],[5,6],[3,4],[5,4],...]
        newConn will be set to None if the sorting algorithm fails.
        """
        barsConn = self.conn
        # We will solve this in an iterative process until the connectivities do not change
        # anymore

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
        self.newConn = newConnFE
        self.newMap = newMap

    def detectTeAngle(self):
        v1 = self.coords[1:-1, :] - self.coords[:-2, :]
        v2 = self.coords[2:, :] - self.coords[1:-1, :]
        v1n = np.linalg.norm(v1, axis=1)
        v2n = np.linalg.norm(v2, axis=1)
        # find the angles between adjacent vectors
        thetas = np.arccos(np.minimum(np.ones(len(self.x) - 2), np.einsum("ij,ij->i", v1, v2) / (v1n[:] * v2n[:])))
        ids = np.argwhere(thetas > 40.0 * np.pi / 180.0)
        if len(ids) == 0:
            print("WARNING, no sharp corner detected")
        elif len(ids) == 1:
            self.lower_te_ind = ids[0][0]
            self.upper_te_ind = ids[0][0]
        elif len(ids) == 2:
            self.lower_te_ind = ids[1][0]
            self.upper_te_ind = ids[0][0]
        else:
            print("WARNING, there are more than 2 sharp corners on slice")
        self.te = 0.5 * (self.coords[self.upper_te_ind, :] + self.coords[self.lower_te_ind, :])

    def detectLE(self):
        te_to_pt_vec = self.coords[:, :] - self.te[:]
        # now compute the distances
        distanmces_to_te = np.linalg.norm(te_to_pt_vec, axis=1)

        # get the max-distance location
        max_dist_ind = np.argmax(distanmces_to_te)

        # roll the values based on this
        for k, v in data.items():
            self.data[k] = np.roll(v, -max_dist_ind)
        coords_array = np.roll(coords_array, -max_dist_ind, axis=0)

        # also adjust the upper and lower TE indices
        upper_te_ind = (upper_te_ind - max_dist_ind) % len(self.data)
        lower_te_ind = (lower_te_ind - max_dist_ind) % len(self.data)

        # save the original LE and TE coordinates for comparison
        orig_le = coords_array[0, :]

        # the LE node and the n-neighboring nodes are interpolated regions with a B-spline
        n_neigh = 10  # number of neighboring nodes to include in both directions
        pts = np.zeros((n_neigh * 2 + 1, 3))

        # center point is always the LE
        pts[n_neigh, :] = coords_array[0, :]

        for jj in range(n_neigh):
            # add the nodes up from the LE
            pts[n_neigh + jj + 1, :] = coords_array[jj + 1, :]
            # add the nodes down from the LE
            pts[n_neigh - jj - 1, :] = coords_array[-1 - jj, :]

        # fit curve to it
        t = np.linspace(0, 1, num=n_neigh * 2 + 1)
        tckx = splrep(t, pts[:, 0])
        tcky = splrep(t, pts[:, 1])
        tckz = splrep(t, pts[:, 2])

        def min_dist_from_te(t):
            vec = np.array([splev(t, tckx) - te[0], splev(t, tcky) - te[1], splev(t, tckz) - te[2]])
            return -np.linalg.norm(vec)

        # can be improved by providing analytic or CS derivatives
        opt_res = minimize(min_dist_from_te, [0], bounds=[[0, 1]], options={"disp": False}, tol=1e-14)

        if not opt_res["success"]:
            print(f"Optimization failed to find the max distance pt at islice {islice}")

        # take the result
        t_opt = opt_res["x"]

        # overwrite the current LE coordinates with the new result
        data["x"][0] = splev(t_opt, tckx)
        data["y"][0] = splev(t_opt, tcky)
        data["z"][0] = splev(t_opt, tckz)
        coords_array[0, :] = np.array([data["x"][0], data["y"][0], data["z"][0]])

        # distance from the old LE to the new LE
        distold = np.linalg.norm(coords_array[0, :] - orig_le)
        # distance from the new LE to upper old LE
        distup = np.linalg.norm(coords_array[0, :] - coords_array[1, :])
        # distance from the old LE to the new LE
        distdown = np.linalg.norm(coords_array[0, :] - coords_array[-1, :])

        if distup >= distdown:
            dist1 = distdown
            dist2 = distold
            interp_ind = -1
        else:
            dist1 = distup
            dist2 = distold
            interp_ind = 1

        dist_tot = dist1 + dist2
        # interpolate all data at the new LE point
        for k, v in data.items():
            # except for coordinates
            if k not in ["x", "y", "z"]:
                old_le_data = v[0]
                interp_pt_data = v[interp_ind]
                # interpolate the new data linearly based on distances.
                data[k][0] = old_le_data * dist1 / dist_tot + interp_pt_data * dist2 / dist_tot

    def plotAirfoil(self):
        fig = go.Figure()
        fig.add_trace(go.Scatter3d(x=self.x, y=self.y, z=self.z, mode="markers", marker=dict(color="black", size=2)))
        if len(self.te) == 3:
            fig.add_trace(
                go.Scatter3d(
                    x=[self.te[0]], y=[self.te[1]], z=[self.te[2]], mode="markers", marker=dict(color="blue", size=4)
                )
            )
        fig.update_layout(
            scene=dict(
                xaxis=dict(color="rgba(0, 0, 0,0)", backgroundcolor="rgba(0, 0, 0,0)"),
                yaxis=dict(color="rgba(0, 0, 0,0)", backgroundcolor="rgba(0, 0, 0,0)"),
                zaxis=dict(color="rgba(0, 0, 0,0)", backgroundcolor="rgba(0, 0, 0,0)"),
            )
        )
        fig.update_scenes(aspectmode="data")
        fig.show()


air = Airfoil(slice_data[0], slice_conn[0])
air.detectTeAngle()
air.plotAirfoil()


# now we are ready to loop over each slice!
for ii, islice in enumerate(range(len(slice_conn))):
    # find the LE x index
    ind_x_min = np.argmin(data["x"])

    # roll all arrays to have LE be the first entry
    for k, v in data.items():
        data[k] = np.roll(v, -ind_x_min)

    # check if the orientation is right, we want the upper skin first (just an arbitrary convention)
    if data["z"][1] - data["z"][0] > 0:
        # we are going up in the LE, this is what we want
        pass
    else:
        for k, v in data.items():
            # flip and roll by 1
            data[k] = np.roll(np.flip(v), 1)

    v1 = coords_array[1:-1, :] - coords_array[:-2, :]
    v2 = coords_array[2:, :] - coords_array[1:-1, :]
    v1n = np.linalg.norm(v1, axis=1)
    v2n = np.linalg.norm(v2, axis=1)
    # find the angles between adjacent vectors
    thetas = np.arccos(np.minimum(np.ones(len(data["x"]) - 2), np.einsum("ij,ij->i", v1, v2) / (v1n[:] * v2n[:])))
    ids = np.argwhere(thetas > 40.0 * np.pi / 180.0)
    if len(ids) == 0:
        print("WARNING, no sharp corner detected")
    elif len(ids) == 1:
        lower_te_ind = ids[0][0] + 1
        upper_te_ind = ids[0][0] + 1
    elif len(ids) == 2:
        lower_te_ind = ids[1][0] + 1
        upper_te_ind = ids[0][0] + 1
    else:
        print("WARNING, there are more than 2 sharp corners on slice")
    # now that we have the TE coordinates, we can find the _true_ LE location,
    # which is the max distance node from the mid-TE point
    te = 0.5 * (coords_array[upper_te_ind, :] + coords_array[lower_te_ind, :])
    # get the vector from the TE to all of the other nodes
    te_to_pt_vec1 = np.zeros_like(coords_array)
    te_to_pt_vec = coords_array[:, :] - te[:]

    # now compute the distances
    distanmces_to_te = np.linalg.norm(te_to_pt_vec, axis=1)

    # get the max-distance location
    max_dist_ind = np.argmax(distanmces_to_te)
    # print(max_dist_ind, len(data["x"]))

    # roll the values based on this
    for k, v in data.items():
        data[k] = np.roll(v, -max_dist_ind)
    coords_array = np.roll(coords_array, -max_dist_ind, axis=0)

    # also adjust the upper and lower TE indices
    upper_te_ind = (upper_te_ind - max_dist_ind) % len(data["x"])
    lower_te_ind = (lower_te_ind - max_dist_ind) % len(data["x"])

    # save the original LE and TE coordinates for comparison
    orig_le = coords_array[0, :]

    # the LE node and the n-neighboring nodes are interpolated regions with a B-spline
    n_neigh = 10  # number of neighboring nodes to include in both directions
    pts = np.zeros((n_neigh * 2 + 1, 3))

    # center point is always the LE
    pts[n_neigh, :] = coords_array[0, :]

    for jj in range(n_neigh):
        # add the nodes up from the LE
        pts[n_neigh + jj + 1, :] = coords_array[jj + 1, :]
        # add the nodes down from the LE
        pts[n_neigh - jj - 1, :] = coords_array[-1 - jj, :]

    # fit curve to it
    t = np.linspace(0, 1, num=n_neigh * 2 + 1)
    tckx = splrep(t, pts[:, 0])
    tcky = splrep(t, pts[:, 1])
    tckz = splrep(t, pts[:, 2])

    def min_dist_from_te(t):
        vec = np.array([splev(t, tckx) - te[0], splev(t, tcky) - te[1], splev(t, tckz) - te[2]])
        return -np.linalg.norm(vec)

    # can be improved by providing analytic or CS derivatives
    opt_res = minimize(min_dist_from_te, [0], bounds=[[0, 1]], options={"disp": False}, tol=1e-14)

    if not opt_res["success"]:
        print(f"Optimization failed to find the max distance pt at islice {islice}")

    # take the result
    t_opt = opt_res["x"]

    # overwrite the current LE coordinates with the new result
    data["x"][0] = splev(t_opt, tckx)
    data["y"][0] = splev(t_opt, tcky)
    data["z"][0] = splev(t_opt, tckz)
    coords_array[0, :] = np.array([data["x"][0], data["y"][0], data["z"][0]])

    # distance from the old LE to the new LE
    distold = np.linalg.norm(coords_array[0, :] - orig_le)
    # distance from the new LE to upper old LE
    distup = np.linalg.norm(coords_array[0, :] - coords_array[1, :])
    # distance from the old LE to the new LE
    distdown = np.linalg.norm(coords_array[0, :] - coords_array[-1, :])

    if distup >= distdown:
        dist1 = distdown
        dist2 = distold
        interp_ind = -1
    else:
        dist1 = distup
        dist2 = distold
        interp_ind = 1

    dist_tot = dist1 + dist2
    # interpolate all data at the new LE point
    for k, v in data.items():
        # except for coordinates
        if k not in ["x", "y", "z"]:
            old_le_data = v[0]
            interp_pt_data = v[interp_ind]
            # interpolate the new data linearly based on distances.
            data[k][0] = old_le_data * dist1 / dist_tot + interp_pt_data * dist2 / dist_tot

    # separate the data for upper and lower parts.
    upper_data = {}
    lower_data = {}
    for k, v in data.items():
        upper_data[k] = v[: upper_te_ind + 1].copy()
        tmp = v[lower_te_ind:].copy()
        # add the LE node to the lower data
        tmp = np.concatenate([tmp, np.array([v[0]])])
        lower_data[k] = tmp
    # now the upper data has all of the nodal data from the LE to the upper TE in an ordered way. the lower data has the same data from the lower TE to the LE.

    # compute the parametric coordinates for both upper and lower data
    for data_dict in [upper_data, lower_data]:
        # vector between nodes on this surface
        vec = np.zeros((len(data_dict["x"]) - 1, 3))
        vec[:, 0] = np.subtract(data_dict["x"][1:], data_dict["x"][:-1])
        vec[:, 1] = np.subtract(data_dict["y"][1:], data_dict["y"][:-1])
        vec[:, 2] = np.subtract(data_dict["z"][1:], data_dict["z"][:-1])
        # distances between all nodes on this surface
        dists = np.linalg.norm(vec, axis=1)
        # initialize the param array to 1 larger to have zero in the beginning
        param = np.zeros(len(dists) + 1)
        # compute the cumulative distance, we put 0 on the first node
        param[1:] = np.cumsum(dists)
        # normalize by total length
        param /= param[-1]
        # save to the original dictionary
        data_dict["param"] = param


def plot_wing():
    fig = go.Figure()
    X = data["x"]
    Y = data["y"]
    Z = data["z"]
    fig.add_trace(go.Scatter3d(x=X, y=Y, z=Z, mode="markers", marker=dict(color="black", size=2)))
    fig.add_trace(go.Scatter3d(x=[te[0]], y=[te[1]], z=[te[2]], mode="markers", marker=dict(color="red", size=3)))

    fig.update_layout(
        scene=dict(
            xaxis=dict(color="rgba(0, 0, 0,0)", backgroundcolor="rgba(0, 0, 0,0)"),
            yaxis=dict(color="rgba(0, 0, 0,0)", backgroundcolor="rgba(0, 0, 0,0)"),
            zaxis=dict(color="rgba(0, 0, 0,0)", backgroundcolor="rgba(0, 0, 0,0)"),
        )
    )
    fig.update_scenes(aspectmode="data")
    fig.show()


plot_wing()
