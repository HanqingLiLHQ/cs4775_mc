# final_graph.py
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os


def generate_final_composite_graph(base_dir, all_cluster_paths):
    composite_index = 1
    for i in range(0, len(all_cluster_paths), 5):
        sub_paths = all_cluster_paths[i : i + 5]
        rows = len(sub_paths)
        fig, axes = plt.subplots(rows, 1, figsize=(8, rows * 4))

        if rows == 1:
            axes = [axes]

        for ax, cluster_graph_path in zip(axes, sub_paths):
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
        final_composite_path = os.path.join(
            base_dir, f"final_composite_{composite_index}.png"
        )
        plt.savefig(final_composite_path, dpi=600)
        plt.close(fig)

        composite_index += 1  # Increment the file index for each composite

    print("All final composite graphs have been generated and saved.\n\n")
