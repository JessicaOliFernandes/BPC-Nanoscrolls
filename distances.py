import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12,           
    "axes.titlesize": 14,      
    "axes.labelsize": 14,      
    "legend.fontsize": 9,      
    "xtick.labelsize": 12,     
    "ytick.labelsize": 12,
})

def plot_nanoscroll_analysis(plane, center, atoms_per_layer, angular_limits, winding_axis, avg_radii, pairwise_distances):
    """
    Generates a plot where each layer is represented by a unique color and marker.
    Distance values (using true nearest-neighbor calculations) are displayed in the legend.
    """
    
    fig, ax = plt.subplots(figsize=(8, 4))
    
    # 1. Plot the original cross-section in light gray
    ax.scatter(plane[:, 0], plane[:, 1], s=1, alpha=0.1, color='gray')
    
    # 2. Plot the REFINED geometric center
    ax.plot(center[0], center[1], color='black', marker='+', markersize=8, mew=2, label='Center (Core)', zorder=5)
    
    # 3. Define Colorblind-Friendly Palette (Okabe-Ito)
    cb_palette = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
    marker_styles = ['o', 's', '^', 'D', 'v', 'p', '*', 'h', 'X', 'd']
    
    core_line_color = '#0072B2' # Okabe-Ito Dark Blue
    gap_line_color = '#D55E00'  # Okabe-Ito Vermillion

    # 4. Plot the classified atoms
    for i, layer in enumerate(atoms_per_layer):
        if not layer: continue
        pts = np.array(layer)
        color = cb_palette[i % len(cb_palette)]
        marker = marker_styles[i % len(marker_styles)]
        label = f'Layer {i+1}' if i < 8 else None
        
        ax.scatter(pts[:, 0], pts[:, 1], s=10, color=color, marker=marker, label=label, zorder=3)

    # 5. Add structural indicators and legends
    dist_legend_handles = []

    if len(avg_radii) > 0:
        core_radius = avg_radii[0]
        base_angle = np.pi / (4)  
        angle_step = 1 * (np.pi / 180)  
        
        # --- CORE RADIUS ---
        x_core = center[0] + core_radius * np.cos(base_angle)
        y_core = center[1] + core_radius * np.sin(base_angle)
        ax.plot([center[0], x_core], [center[1], y_core], color=core_line_color, linestyle=':', linewidth=1.5, zorder=4)
        
        dist_legend_handles.append(Line2D([0], [0], color=core_line_color, linestyle=':', linewidth=1.5, 
                                          label=rf"$r_{{core}} = {core_radius:.2f}$ Å"))

        # --- PAIRWISE INTER-LAYER DISTANCES (Nearest Neighbor) ---
        if len(pairwise_distances) > 0:
            for j in range(len(pairwise_distances)): 
                if j >= 6: break 
                
                r_start = avg_radii[j]
                r_end = avg_radii[j+1]
                dist_val = pairwise_distances[j]
                
                current_angle = base_angle + ((j+1) * angle_step)
                
                x_start = center[0] + r_start * np.cos(current_angle)
                y_start = center[1] + r_start * np.sin(current_angle)
                x_end = center[0] + r_end * np.cos(current_angle)
                y_end = center[1] + r_end * np.sin(current_angle)
                
                ax.plot([x_start, x_end], [y_start, y_end], color=gap_line_color, linestyle=':', linewidth=1.5, zorder=4)
                
                dist_legend_handles.append(Line2D([0], [0], color=gap_line_color, linestyle=':', linewidth=1.5, 
                                                  label=rf"$d_{{{j+1}-{j+2}}} = {dist_val:.2f}$ Å"))
            
            mean_overall_distance = np.mean(pairwise_distances)
            dist_legend_handles.append(Line2D([0], [0], color='none', 
                                              label=rf"Mean $d = {mean_overall_distance:.2f}$ Å"))

    # 6. Formatting the plot
    if winding_axis.lower() == 'x':
        xlabel, ylabel = "Y", "Z"
    elif winding_axis.lower() == 'y':
        xlabel, ylabel = "X", "Z"
    else:
        xlabel, ylabel = "X", "Y"

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel(f"Coordinate {xlabel} (Å)")
    ax.set_ylabel(f"Coordinate {ylabel} (Å)")
    ax.set_title(f"Nanoscroll Layer Analysis \n Repetition 4 - Zigzag type", pad=10, fontsize=12)

    if len(atoms_per_layer) > 0:
        main_legend = ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0), 
                                frameon=True, edgecolor='black', title="Structure",
                                fancybox=False, fontsize=9, markerscale=1.2)
        ax.add_artist(main_legend) 
    
    if dist_legend_handles:
        ax.legend(handles=dist_legend_handles, loc='lower left', bbox_to_anchor=(1.05, 0.0), 
                  frameon=True, edgecolor='black', title="Distances",
                  fancybox=False, fontsize=9)
    
    plt.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.3)
    
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.show()

def _find_layers(plane, center, n_slices, gap_threshold):
    """Helper function to perform layer detection."""
    dx = plane[:, 0] - center[0]
    dy = plane[:, 1] - center[1]
    
    radii = np.sqrt(dx**2 + dy**2)
    angles = np.arctan2(dy, dx)

    angular_limits = np.linspace(-np.pi, np.pi, n_slices + 1)
    atoms_per_layer = [] 

    for i in range(n_slices):
        mask = (angles >= angular_limits[i]) & (angles < angular_limits[i+1])
        slice_indices = np.where(mask)[0]
        
        if len(slice_indices) == 0: continue
            
        sort_order = np.argsort(radii[slice_indices])
        sorted_indices = slice_indices[sort_order]
        
        layer_idx = 0
        r_previous = radii[sorted_indices[0]]
        
        for idx in sorted_indices:
            r_current = radii[idx]
            
            if (r_current - r_previous) > gap_threshold:
                layer_idx += 1
            
            while len(atoms_per_layer) <= layer_idx:
                atoms_per_layer.append([])
            
            atoms_per_layer[layer_idx].append(plane[idx])
            r_previous = r_current
            
    return atoms_per_layer, angular_limits

def analyze_nanoscroll_radii(file_coords, winding_axis='x', n_slices=12, gap_threshold=1.8):
    """Calculates true nearest-neighbor layer distances robust to partial layers."""
    try:
        coords = np.loadtxt(file_coords, skiprows=2, usecols=(1, 2, 3))
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    if winding_axis.lower() == 'x':
        plane = coords[:, [1, 2]] 
    elif winding_axis.lower() == 'y':
        plane = coords[:, [0, 2]] 
    else:
        plane = coords[:, [0, 1]] 

    initial_center = np.mean(plane, axis=0)
    initial_layers, _ = _find_layers(plane, initial_center, n_slices, gap_threshold)
    
    if len(initial_layers) > 0 and len(initial_layers[0]) > 0:
        refined_center = np.mean(np.array(initial_layers[0]), axis=0)
    else:
        refined_center = initial_center 
        
    atoms_per_layer, angular_limits = _find_layers(plane, refined_center, n_slices, gap_threshold)

    avg_radii = []
    pairwise_distances = []
    
    output_text = f"=== Results for Nanoscroll ({winding_axis.upper()}-axis) ===\n"
    output_text += f"Initial 2D Center Guess: X={initial_center[0]:.3f} Å, Y={initial_center[1]:.3f} Å\n"
    output_text += f"Refined 2D Center (Core): X={refined_center[0]:.3f} Å, Y={refined_center[1]:.3f} Å\n\n"
    
    if len(atoms_per_layer) > 0:
        for layer in atoms_per_layer:
            if not layer: continue
            layer_arr = np.array(layer)
            r_layer = np.sqrt((layer_arr[:,0]-refined_center[0])**2 + (layer_arr[:,1]-refined_center[1])**2)
            avg_radii.append(np.mean(r_layer))
        
        output_text += f"Average CORE Radius (Layer 1): {avg_radii[0]:.3f} Å\n\n"
        
        if len(atoms_per_layer) > 1:
            output_text += "True Nearest-Neighbor Inter-Layer Distances:\n"
            
            for j in range(len(atoms_per_layer) - 1):
                pts1 = np.array(atoms_per_layer[j])
                pts2 = np.array(atoms_per_layer[j+1])
                
                if len(pts1) == 0 or len(pts2) == 0:
                    continue
                
                diff = pts1[:, np.newaxis, :] - pts2[np.newaxis, :, :]
                dist_matrix = np.sqrt(np.sum(diff**2, axis=-1))
                
                min_dists_1_to_2 = np.min(dist_matrix, axis=1) 
                min_dists_2_to_1 = np.min(dist_matrix, axis=0) 
                
                all_min_dists = np.concatenate([min_dists_1_to_2, min_dists_2_to_1])
                base_gap = np.percentile(all_min_dists, 10)
                
                valid_dists = all_min_dists[all_min_dists < base_gap * 1.4] 
                
                mean_d = np.mean(valid_dists)
                std_d = np.std(valid_dists)
                
                pairwise_distances.append(mean_d)
                
                output_text += f"  Layer {j+1} -> Layer {j+2}: {mean_d:.3f} Å (± {std_d:.3f})\n"
            
            if pairwise_distances:
                overall_mean = np.mean(pairwise_distances)
                output_text += f"\nOverall Average Inter-Layer Distance: {overall_mean:.3f} Å\n"
            
    output_text += "==========================================\n"
    
    print(output_text) 
    
    output_filename = "distancias_nanoscroll.txt"
    try:
        with open(output_filename, 'w', encoding='utf-8') as f:
            f.write(output_text)
        print(f"-> Data successfully saved to '{output_filename}'\n")
    except Exception as e:
        print(f"-> Could not save to file: {e}\n")
    
    plot_nanoscroll_analysis(plane, refined_center, atoms_per_layer, angular_limits, winding_axis, avg_radii, pairwise_distances)

if __name__ == "__main__":
    analyze_nanoscroll_radii("geo.final.xyz", n_slices=12)
