import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from Bio import motifs
import os


class MotifGraphGenerator:
    def __init__(self, dataset_cc, base_dir):
        self.dataset_cc = dataset_cc
        self.base_dir = base_dir
        os.makedirs(self.base_dir, exist_ok=True)

    def save_motif_logo(self, motif, filename):
        motif.weblogo(filename)

    def generate_cluster_graphs(self):
        all_cluster_paths = []

        for i, motif_dict in enumerate(self.dataset_cc, start=1):
            cluster_dir = os.path.join(self.base_dir, f"cluster_{i}")
            os.makedirs(cluster_dir, exist_ok=True)
            image_paths = []

            for motif_name, motif_obj in motif_dict.items():
                logo_filename = f"{motif_name}.png"
                logo_path = os.path.join(cluster_dir, logo_filename)
                image_paths.append(logo_path)
                self.save_motif_logo(motif_obj, logo_path)

            num_motifs = len(motif_dict)
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

            plt.tight_layout()
            cluster_graph_path = os.path.join(cluster_dir, f"cluster_{i}.png")
            plt.savefig(cluster_graph_path, dpi=600)
            plt.close(fig)
            all_cluster_paths.append(cluster_graph_path)

        return all_cluster_paths

    def generate_final_composite_graph(self, all_cluster_paths):
        total_clusters = len(all_cluster_paths)
        final_fig, final_axes = plt.subplots(
            total_clusters, 1, figsize=(8, total_clusters * 4)
        )

        if total_clusters == 1:
            final_axes = [final_axes]

        for ax, cluster_graph_path in zip(final_axes, all_cluster_paths):
            try:
                cluster_img = mpimg.imread(cluster_graph_path)
                ax.imshow(cluster_img)
                ax.axis("off")
                cluster_name = os.path.basename(cluster_graph_path).replace(".png", "")
                ax.set_title(cluster_name)
            except FileNotFoundError:
                print(f"File {cluster_graph_path} not found. Skipping.")
            except SyntaxError:
                print(f"File {cluster_graph_path} is not a valid PNG. Skipping.")

        plt.tight_layout()
        final_composite_path = os.path.join(self.base_dir, "final_composite.png")
        plt.savefig(final_composite_path, dpi=600)
        plt.close(final_fig)

        print("Final composite graph has been generated and saved.")
