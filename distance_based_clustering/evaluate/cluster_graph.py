# graph.py
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os


class ClusterGraphGenerator:
    def __init__(self, clusters, motifs_dir, output_dir):
        self.clusters = (
            clusters  # List of lists, where each inner list contains motif names
        )
        self.motifs_dir = motifs_dir  # Directory containing individual motif images
        self.output_dir = output_dir  # Output directory for cluster graphs
        os.makedirs(self.output_dir, exist_ok=True)

    def generate_cluster_graphs(self):
        all_cluster_paths = []

        for i, cluster in enumerate(self.clusters, start=1):
            cluster_dir = os.path.join(self.output_dir, f"cluster_{i}")
            os.makedirs(cluster_dir, exist_ok=True)
            image_paths = [
                os.path.join(self.motifs_dir, f"{motif_name}.png")
                for motif_name in cluster
            ]

            num_motifs = len(cluster)

            if num_motifs <= 5:
                fig, axes = plt.subplots(1, num_motifs, figsize=(num_motifs * 4, 4))
                if num_motifs == 1:
                    axes = [axes]
                for ax, logo_path in zip(axes, image_paths):
                    try:
                        img = mpimg.imread(logo_path)
                        ax.imshow(img)
                        ax.axis("off")
                        motif_name = os.path.basename(logo_path).replace(".png", "")
                        ax.set_title(motif_name)
                    except FileNotFoundError:
                        print(f"File {logo_path} not found. Skipping.")
                    except SyntaxError:
                        print(f"File {logo_path} is not a valid PNG. Skipping.")
            else:
                rows = (
                    num_motifs + 4
                ) // 5  # Calculate how many rows are needed for 5 logos per row
                cols = min(num_motifs, 5)  # Use up to 5 columns
                fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 4))
                axes = axes.flatten() if rows > 1 else [axes]  # Flatten if necessary
                for ax in axes:
                    ax.axis("off")  # Hide all axes by default
                for ax, logo_path in zip(axes, image_paths):
                    try:
                        img = mpimg.imread(logo_path)
                        ax.imshow(img)
                        ax.axis("off")  # Only show axes for images
                        motif_name = os.path.basename(logo_path).replace(".png", "")
                        ax.set_title(motif_name)
                    except FileNotFoundError:
                        print(f"File {logo_path} not found. Skipping.")
                    except SyntaxError:
                        print(f"File {logo_path} is not a valid PNG. Skipping.")

            plt.tight_layout()
            cluster_graph_path = os.path.join(cluster_dir, f"cluster_{i}.png")
            plt.savefig(cluster_graph_path, dpi=600)
            plt.close(fig)
            all_cluster_paths.append(cluster_graph_path)

        print("All cluster graphs have been generated and saved. \n\n")
        return all_cluster_paths
