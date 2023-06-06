import numpy as np
import warnings
import matplotlib.pyplot as plt


def findElem(elems, index, val):
    for i in range(len(elems)):
        if elems[i][index] == val:
            return i

    return None


def FEsort(barsConn):
    """
    This function can be used to sort connectivities coming from the CGNS file

    barsConn should be a list of integers [[2,3],[5,6],[3,4],[5,4],...]

    newConn will be set to None if the sorting algorithm fails.
    """

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
    return newConnFE, newMap


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
    slice_y = np.zeros(nmax)

    # beginning line index for the current slice
    slice_begin = 2

    # loop over slices
    for islice in range(nmax):
        slice_header = lines[slice_begin].replace('"', "").replace("(", "").replace(")", "").split()
        y_loc = float(slice_header[-1])

        slice_y[islice] = y_loc

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
        # print(tmp["CoordinateX"])
        slice_data.append(tmp)

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

        # use FESort to fix the conn
        newConn, newMap = FEsort(conn.tolist())

        # the line segments in newConn is likely unordered, we need to re-order them
        x = tmp["CoordinateX"]
        y = tmp["CoordinateY"]
        z = tmp["CoordinateZ"]

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

        # save this in the large list
        slice_conn.append(np.array(final_conn))

        # increment the beginning line index for the next slice
        # + 1 is for the Datapacking=point line, + 2 is for the 2 headers
        slice_begin += 1 + n_node + n_elem + 2

    data_dict = {
        "slice_data": slice_data,
        "slice_conn": slice_conn,
        "slice_y": slice_y,
    }

    return data_dict


def readLift(filename):
    # get all lines
    with open(filename, "r") as f:
        lines = f.readlines()

    # first line is the title
    # second line has the variable list
    var_list = (
        lines[1].replace('" "', "$").replace('"', "").replace("=", "$").replace(" ", "").replace("\n", "").split("$")
    )
    var_list.remove("Variables")
    print("Found variables:", var_list)

    # number of variables
    n_var = len(var_list)

    # get the number of stations
    n_pt = int(lines[3].split()[-1])

    # read the data in one long array
    tmp = np.genfromtxt(
        filename,
        dtype={
            "names": ["data"],
            "formats": ["f4"],
        },
        skip_header=5,
        max_rows=n_pt * n_var,
    )

    # loop over each data and extract values to a dict
    data_dict = {}
    for ii, var in enumerate(var_list):
        data_dict[var] = np.array(tmp[ii * n_pt : (ii + 1) * n_pt])

    return data_dict
