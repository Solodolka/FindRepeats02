from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
from Bio.Align import substitution_matrices
import numpy as np
import re

genomes_dir = "genomes"
repeats_dir = "repeats"
results_dir = "results"
species_dir = "beet"
genomes_amount = 7

def create_custom_blastn_matrix(match_score=2, mismatch_score=-1, n_vs_n_score=2, n_vs_acgt_score=0.5):
    """
    Создает кастомную матрицу скоринга, аналогичную BLASTN,
    с возможностью корректировки значений для 'N', возвращая объект Array.

    Args:
        match_score (int): Оценка за совпадение A-A, T-T, G-G, C-C.
        mismatch_score (int): Оценка за несовпадение A-G, T-A и т.д.
        n_vs_n_score (int): Оценка за совпадение N-N.
        n_vs_acgt_score (int): Оценка за несовпадение N-A, N-T, N-G, N-C и наоборот.

    Returns:
        substitution_matrices.Array: Кастомная матрица скоринга для PairwiseAligner.
    """
    alphabet = "ACGTNRYMKSW"
    letter_to_index = {letter: i for i, letter in enumerate(alphabet)}
    matrix_data = np.zeros((len(alphabet), len(alphabet)), dtype=float)

    for char1 in alphabet:
        for char2 in alphabet:
            i = letter_to_index[char1]
            j = letter_to_index[char2]
            if char1 == char2:
                if char1 == 'N':
                    matrix_data[i, j] = n_vs_n_score
                else:
                    matrix_data[i, j] = match_score
            else:
                if char1 == 'N' or char2 == 'N':
                    matrix_data[i, j] = n_vs_acgt_score
                else:
                    matrix_data[i, j] = mismatch_score

    return substitution_matrices.Array(alphabet = alphabet, data = matrix_data)


def find_and_write_fragment_positions(gene_record, fragment_records):
    """ Функция ищет фрагменты в гене и записывает результаты в соответствующие файлы. :param gene_record: запись гена :param fragment_records: список записей фрагментов """
    gene_id = gene_record.id
    gene_seq = gene_record.seq.upper()
    seq1 = gene_seq

    for fragment_record in fragment_records:
        fragment_id = fragment_record.id
        fragment_seq = fragment_record.seq.upper()
        seq2 = fragment_seq

        # custom_matrix = create_custom_blastn_matrix()
        custom_matrix = substitution_matrices.load("NUC.4.4")

        aligner = Align.PairwiseAligner(open_gap_score=-11.0,
                                        extend_gap_score=-2.0,
                                        mode='local',
                                        substitution_matrix = custom_matrix)
        print(aligner.algorithm)

        alignments = aligner.align(seq1, seq2)
        # min_score = len(seq2) * 1.7
        # filtered_alignments_iterator = filter(lambda alignment: alignment.score >= min_score, alignments)
        unic_alignments = {}
        for al in alignments:
            key = str (al.coordinates[0][0])
            if key in unic_alignments:
                if al.score > unic_alignments[key].score:
                    unic_alignments[key] = al
            else:
                unic_alignments[key] = al

        print('Alignmets number: ', len(alignments))
        if len(unic_alignments) > 0:
            # Создаем файл для вывода результатов
            fragment_filename = re.sub(r'\W+', '', fragment_id)
            fragment_filename = f"{results_dir}/{species_dir}/{fragment_filename}.fasta"
            with open(fragment_filename, 'a') as output_file:
                for al in unic_alignments.values():
                    print('Normalized score: ', al.score / len(seq2))
                    # print('Alignmet coordinates: ')
                    # print(al.coordinates)
                    # print('Alignmet substitutions: ')
                    # print(al.substitutions)
                    print(f"Repeat {fragment_id} found in chromosome {gene_id} on position {al.coordinates}")

                    # Формируем строку заголовка
                    header = f">{gene_id}\tPosition: {al.coordinates}\tScore: {al.score / len(seq2)}"
                    header = re.sub(r'\n', '', header)

                    # Определяем контекст вокруг фрагмента
                    start_pos = max(0, al.coordinates[0][0] - 31)
                    arr_length = len(al.coordinates[0])
                    end_pos = min(len(gene_seq), al.coordinates[0][arr_length-1] + 29)

                    # Формируем строку контекста
                    context = gene_seq[start_pos:end_pos]

                    # Записываем результат в файл
                    output_file.write(header + '\n' + str(context) + '\n')
        else:
            print(f"Repeat {fragment_id} not found in chromosome {gene_id}")


# Чтение файла с фрагментами
fragments_path = f"{repeats_dir}/{species_dir}/Repeats.fasta"
fragment_records = list(SeqIO.parse(fragments_path, "fasta"))
print('Repeats are loaded')

for file_num in range(1, genomes_amount+1):
    # Чтение файла с геномом
    gene_path = f"{genomes_dir}/{species_dir}/Gene{file_num}.fasta"
    try:
        gene_record = next(SeqIO.parse(gene_path, "fasta"))
    except:
        raise ValueError(f"Genom reading erroe in {gene_path}")
    print(f"Chromosome {gene_path} is loaded")

    # Поиск и запись результатов
    find_and_write_fragment_positions(gene_record, fragment_records)
