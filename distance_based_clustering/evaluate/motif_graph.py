# motif_graph.py
from Bio import motifs
import os


class MotifGraphGenerator:
    def __init__(self, dataset, output_dir):
        self.dataset = (
            dataset  # This should be a dictionary of {motif_name: Motif_object}
        )
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def save_motif_logo(self, motif, filename):
        motif.weblogo(filename)

    def generate_motif_graphs(self):
        for motif_name, motif_obj in self.dataset.items():
            logo_filename = f"{motif_name}.png"
            logo_path = os.path.join(self.output_dir, logo_filename)
            self.save_motif_logo(motif_obj, logo_path)

        print("All motif graphs have been generated and saved. \n\n")
