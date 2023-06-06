import numpy as np
import warnings
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from parseSlice import *
import plotly.graph_objects as go

slice_dict = readSlices("adflow/slice.dat", nmax=-1)

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
    te_to_pt_vec = np.zeros_like(coords_array)

    for idim in range(3):
        te_to_pt_vec[:, idim] = np.subtract(coords_array[:, idim], mid_te[idim])

    # now compute the distances
    distanmces_to_te = np.linalg.norm(te_to_pt_vec, axis=1)

    # get the max-distance location
    max_dist_ind = np.argmax(distanmces_to_te)
    # print(max_dist_ind, len(data["x"]))

    # roll the values based on this
    for k, v in data.items():
        data[k] = np.roll(v, -max_dist_ind)

    # also adjust the upper and lower TE indices
    upper_te_ind = (upper_te_ind - max_dist_ind) % len(data["x"])
    lower_te_ind = (lower_te_ind - max_dist_ind) % len(data["x"])

    # save the original LE and TE coordinates for comparison
    orig_le_z = data["z"][0]
    orig_le_x = data["x"][0]

    # the LE node and the n-neighboring nodes are interpolated regions with a p_poly order polynomial
    # TODO using more neighbors will probably be more accurate.
    n_neigh = 1  # number of neighboring nodes to include in both directions
    p_poly = 2  # order of the polynomial fit
    pts = np.zeros((n_neigh * 2 + 1, 2))

    # center point is always the LE
    pts[n_neigh, 0] = data["z"][0]
    pts[n_neigh, 1] = data["x"][0]

    for jj in range(n_neigh):
        # add the nodes up from the LE
        pts[n_neigh + jj + 1, 0] = data["z"][jj + 1]
        pts[n_neigh + jj + 1, 1] = data["x"][jj + 1]

        # add the nodes down from the LE
        pts[n_neigh - jj - 1, 0] = data["z"][-1 - jj]
        pts[n_neigh - jj - 1, 1] = data["x"][-1 - jj]

    # fit a 2nd order polynomial
    coefs = np.polyfit(pts[:, 0], pts[:, 1], p_poly)

    def my_poly(z_loc):
        val = 0.0
        for kk in range(p_poly):
            val += coefs[kk] * (z_loc ** (p_poly - kk))
        val += coefs[-1]
        return val

    def min_dist_from_te(z_loc):
        x_loc = my_poly(z_loc)
        vec = np.array([z_loc - mid_te[2], x_loc - mid_te[0]])
        return -np.linalg.norm(vec)

    # bounds are always the first neighbors. we dont want to search any further than them
    bounds = [[pts[n_neigh - 1, 0], pts[n_neigh + 1, 0]]]

    # TODO optimization can be improved by providing analytic or CS derivatives
    opt_res = minimize(min_dist_from_te, pts[n_neigh, 0], bounds=bounds, options={"disp": False}, tol=1e-14)

    if not opt_res["success"]:
        print(f"Optimization failed to find the max distance pt at islice {islice}")

    # take the result
    z_opt = opt_res["x"]
    x_opt = my_poly(z_opt)

    # overwrite the current LE coordinates with the new result
    data["z"][0] = z_opt
    data["x"][0] = x_opt

    # check if the new point is above or below the old one
    if z_opt > orig_le_z:
        # new point is above, so we use the old LE and the first point above it to interpolate
        interp_ind = 1
    else:
        # new point is below the old LE, we use the old LE and the first point below to interpolate
        interp_ind = -1

    # compute the distances from the new LE to the old LE and the other interp point

    # distance from the old LE to the new LE
    dist1 = np.linalg.norm(
        np.array(
            [
                x_opt - orig_le_x,
                z_opt - orig_le_z,
            ]
        )
    )

    # distance from the new LE to the interp point
    dist2 = np.linalg.norm(
        np.array(
            [
                x_opt - data["x"][interp_ind],
                z_opt - data["z"][interp_ind],
            ]
        )
    )

    # total distance
    dist_tot = dist1 + dist2

    # interpolate all data at the new LE point
    for k, v in data.items():
        # except for coordinates
        if k not in ["x", "y", "z"]:
            old_le_data = v[0]
            interp_pt_data = v[interp_ind]
            # interpolate the new data linearly based on distances. crude but should work
            # another option would be just interpolating by z-coordinate
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


class Airfoil:
    def __init__(self):
        pass


def plot_wing():
    fig = go.Figure()
    X = data["x"]
    Y = data["y"]
    Z = data["z"]
    fig.add_trace(go.Scatter3d(x=X, y=Y, z=Z, mode="markers", marker=dict(color="black", size=2)))
    fig.add_trace(
        go.Scatter3d(x=[mid_te[0]], y=[mid_te[1]], z=[mid_te[2]], mode="markers", marker=dict(color="red", size=3))
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


plot_wing()
