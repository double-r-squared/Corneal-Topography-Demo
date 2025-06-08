from mpl_toolkits.mplot3d import Axes3D  # registers 3D projection
import plotly.graph_objects as go


# Flatten all (y, x) points
all_yx = [pt for angle_list in cartesian_coords for pt in angle_list]

def analyze_gaps_per_ray(polar_data):
    smallest_gaps = []
    largest_gaps = []
    median_gaps = []

    for radii in polar_data:
        # Calculate gaps for this ray
        diffs = [radii[i+1] - radii[i] for i in range(len(radii) - 1)]
        
        smallest_gaps.append(min(diffs))
        largest_gaps.append(max(diffs))
        median_gaps.append(np.median(diffs))

    # Now find median of these arrays across all rays
    median_of_smallest = np.median(smallest_gaps)
    median_of_largest = np.median(largest_gaps)
    median_of_medians = np.median(median_gaps)

    print("Across all rays:")
    print(f"Median of smallest gaps: {median_of_smallest}")
    print(f"Median of largest gaps: {median_of_largest}")
    print(f"Median of median gaps: {median_of_medians}")

    return smallest_gaps, largest_gaps, median_gaps

# Example usage:
smallest_gaps, largest_gaps, median_gaps = analyze_gaps_per_ray(averaged_polar_data)

def index_rings_for_ray(radii, min_gap, max_gap):
    indices = [0]  # ring index for the first radius is 0
    for i in range(1, len(radii)):
        gap = radii[i] - radii[i-1]
        if min_gap <= gap <= max_gap:
            # normal ring increment by 1
            indices.append(indices[-1] + 1)
        else:
            # gap too small or too large → jump by 2 (or you can decide differently)
            indices.append(indices[-1] + 2)
    return indices

def index_rings_all_rays(all_rays, min_gap, max_gap):
    all_indices = []
    for ray in all_rays:
        ring_indices = index_rings_for_ray(ray, min_gap, max_gap)
        all_indices.append(ring_indices)
    return all_indices


min_gap = 4.0
max_gap = 8.0

ring_indices_all = index_rings_all_rays(averaged_polar_data, min_gap, max_gap)


# Print ring indices for first ray as example:
print(ring_indices_all)

def assign_indices_to_coords(coords, indices):
    if len(coords) != len(indices):
        raise ValueError("coords and indices must be the same length")
    return [(y, x, idx) for (y, x), idx in zip(coords, indices)]

def assign_indices_to_all_coords(all_coords, all_indices):
    all_indexed_coords = []
    for coords, indices in zip(all_coords, all_indices):
        all_indexed_coords.append(assign_indices_to_coords(coords, indices))
    return all_indexed_coords

indexed_coords = assign_indices_to_all_coords(cartesian_coords, ring_indices_all)
def filter_indexed_coords_by_max_ring(indexed_coords, max_index=32):
    filtered = []
    for ray in indexed_coords:
        filtered_ray = [(y, x, idx) for y, x, idx in ray if idx <= max_index]
        filtered.append(filtered_ray)
    return filtered

indexed_coords = filter_indexed_coords_by_max_ring(indexed_coords, max_index=32)


def plot_indexed_coords_3d_plotly(rays_indexed_coords):
    flat_coords = [pt for ray in rays_indexed_coords for pt in ray if len(pt) == 3]

    ys = [y for y, x, idx in flat_coords]
    xs = [x for y, x, idx in flat_coords]
    zs = [idx for y, x, idx in flat_coords]

    fig = go.Figure(data=[go.Scatter3d(
        x=xs, y=ys, z=zs,
        mode='markers',
        marker=dict(
            size=4,
            color=zs,
            colorscale='Viridis',
            opacity=0.8,
            colorbar=dict(title='Ring Index')
        )
    )])

    fig.update_layout(
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Ring Index'
        ),
        title='3D Scatter Plot of Indexed Placido Ring Coordinates',
        margin=dict(l=0, r=0, b=0, t=30)
    )
    fig.show()

plot_indexed_coords_3d_plotly(indexed_coords)