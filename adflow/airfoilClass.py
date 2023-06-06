import numpy as np
import warnings
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from parseSlice import *
import plotly.graph_objects as go
from scipy.interpolate import splrep, splev


class Airfoil:
    """
    Create the Airfoil object.

    Parameters
    ----------
    data : dict
        dict output by the cut and reordered.
    """

    def __init__(self, data):
        self.data = data
        self.x = data["x"]
        self.y = data["y"]
        self.z = data["z"]
        self.u = data["u"]
        self.v = data["v"]
        self.w = data["w"]
        self.cp = data["cp"]

    def detectTeAngle(self):
        pass


slice_dict = readSlices("input_files/slice.dat", nmax=-1)

slice_data = slice_dict["slice_data"]
slice_conn = slice_dict["slice_conn"]
slice_y = slice_dict["slice_y"]
# check how many slices we will actually output:
n_slice = len(slice_conn)

# now we are ready to loop over each slice!
for ii, islice in enumerate(range(len(slice_conn))):
    # ii is the index of the slice we are writing, islice is the index of the slice in the original dataset
    y_loc = slice_y[islice]

    conn = slice_conn[islice]
    # get the node list
    node_list = conn[:, 0].tolist()

    data = {
        "x": slice_data[islice]["CoordinateX"][node_list],
        "y": slice_data[islice]["CoordinateY"][node_list],
        "z": slice_data[islice]["CoordinateZ"][node_list],
        "u": slice_data[islice]["VelocityX"][node_list],
        "v": slice_data[islice]["VelocityY"][node_list],
        "w": slice_data[islice]["VelocityZ"][node_list],
        "cp": slice_data[islice]["CoefPressure"][node_list],
    }

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

    # loop over elements and find the 2 bends that are larger than X degrees
    upper_te_found = False
    lower_te_found = False

    # put all of the coordinates in one n-node by 3 matrix for easier computation
    coords_array = np.zeros((len(data["x"]), 3))
    coords_array[:, 0] = data["x"].copy()
    coords_array[:, 1] = data["y"].copy()
    coords_array[:, 2] = data["z"].copy()

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
