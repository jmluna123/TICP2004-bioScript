from Bio import AlignIO
from ete3 import Tree
import numpy as np
import sys
import os

extensions = [".fasta", ".fna", ".ffn", ".faa", ".frn"]


def print_help():
    print(
        """>>>>UPGMA<<<<
usage: python UPGMA.py <route>"""
    )
    exit()


def print_error(message):
    print("\033[91m" + message + "\033[0m")
    print_help()


def verify_path(path):
    filename = os.path.splitext(os.path.basename(path))

    if filename[1] not in extensions:
        print_error("ERROR: file should be .fasta format")
    if not os.path.isfile(path):
        print_error("ERROR: please enter a valid file path")

    return filename


def count_diff_sequences(seq1, seq2):
    diffs = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            diffs += 1
    return diffs


def get_distance_matrix(sequences):
    tags = list(sequences)
    n_sequences = len(sequences)
    distance_matrix = np.zeros((n_sequences, n_sequences))

    for i in range(n_sequences):
        for j in range(i + 1, n_sequences):
            distance_matrix[i, j] = distance_matrix[j, i] = count_diff_sequences(
                sequences[tags[i]], sequences[tags[j]]
            )

    return distance_matrix


class UPGMA:
    def __init__(self, sequences, distance_matrix):
        self.sequences = sequences
        self.tmp_sequences = sequences
        self.distance_matrix = distance_matrix
        self.num_sequences = len(sequences)
        self.tree = []

    def run(self):
        while len(self.distance_matrix) > 1:
            i, j = self.find_closest_sequences()
            self.update_distance_matrix(i, j)
            self.update_tree(i, j)

    def find_closest_sequences(self):
        min_distance = np.inf
        min_i, min_j = -1, -1

        for i in range(len(self.distance_matrix)):
            for j in range(i + 1, len(self.distance_matrix[i])):
                if self.distance_matrix[i][j] < min_distance:
                    min_distance = self.distance_matrix[i][j]
                    min_i, min_j = i, j

        return min_i, min_j

    def update_distance_matrix(self, i, j):
        new_matrix = []
        new_row = []
        new_row.append(0.0)
        for k in range(len(self.distance_matrix)):
            if k != i and k != j:
                new_distance = (
                    self.distance_matrix[i][k] + self.distance_matrix[j][k]
                ) / 2
                new_row.append(new_distance)
        new_matrix.append(new_row)

        for x in range(1, len(self.distance_matrix)):
            if x != i and x != j:
                new_row = []
                for y in range(i + 1, len(self.distance_matrix[i])):
                    if y != i and y != j:
                        new_row.append(self.distance_matrix[x][y])
                    else:
                        new_row.append(0.0)
                new_matrix.append(new_row)
        self.distance_matrix = new_matrix

    def update_tree(self, i, j):
        seq_i = self.tmp_sequences[i]
        seq_j = self.tmp_sequences[j]
        node = f"({seq_i},{seq_j})"

        tmp = []
        tmp.append(node)
        for k in range(len(self.tmp_sequences)):
            if k != i and k != j:
                tmp.append(self.tmp_sequences[k])
        self.tmp_sequences = tmp

        self.tree.append(node)


if __name__ == "__main__":
    args = sys.argv[1:]
    n_args = len(args)

    route = ""

    if n_args == 0:
        print_error("ERROR: There is no route specified")

    route = args[0]
    filename = verify_path(route)

    # reed fasta
    sequences = {}
    for record in AlignIO.read(route, "fasta"):
        sequences[record.id] = str(record.seq)

    distance_matrix = get_distance_matrix(sequences)

    print(">>>> Distance Matrix <<<<")
    print(distance_matrix)

    upgma = UPGMA(list(sequences), distance_matrix)
    upgma.run()

    print("Final Tree:", upgma.tree[-1])
    t = Tree(upgma.tree[-1] + ";")
    print(t)
    t.show()
