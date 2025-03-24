from Bio import SeqIO
import re

genomes_dir = "genomes"
repeats_dir = "repeats"
results_dir = "results"
species_dir = "test"
genomes_amount = 2

def find_and_write_fragment_positions(gene_record, fragment_records):
    """ Функция ищет фрагменты в гене и записывает результаты в соответствующие файлы. :param gene_record: запись гена :param fragment_records: список записей фрагментов """
    gene_id = gene_record.id
    gene_seq = gene_record.seq.upper()

    for fragment_record in fragment_records:
        fragment_id = fragment_record.id
        fragment_seq = fragment_record.seq.upper()

        positions = []
        start_index = 0
        while True:
            position = gene_seq.find(fragment_seq, start_index)
            if position == -1:
                break
            positions.append(position + 1)  # Положение в гене начинается с 1
            start_index = position + 1

        if positions:
            print(f"Фрагмент {fragment_id} найден в гене {gene_id} на позициях {positions}")

            # Создаем файл для вывода результатов
            fragment_filename = re.sub(r'\W+', '', fragment_id)
            fragment_filename = f"{results_dir}/{species_dir}/{fragment_filename}.fasta"
            with open(fragment_filename, 'a') as output_file:
                for pos in positions:
                    # Формируем строку заголовка
                    header = f">{gene_id}\tPosition: {pos}"

                    # Определяем контекст вокруг фрагмента
                    start_pos = max(0, pos - 31)
                    end_pos = min(len(gene_seq), pos + len(fragment_seq) + 29)

                    # Формируем строку контекста
                    context = gene_seq[start_pos:pos - 1]
                    context += fragment_seq
                    context += gene_seq[pos + len(fragment_seq) - 1:end_pos]

                    # Записываем результат в файл
                    output_file.write(header + '\n' + str(context) + '\n')
        else:
            print(f"Фрагмент {fragment_id} не найден в гене {gene_id}")


# Чтение файла с фрагментами
fragments_path = f"{repeats_dir}/{species_dir}/Repeats.fasta"
fragment_records = list(SeqIO.parse(fragments_path, "fasta"))
print('Repeats are loaded')

for file_num in range(1, genomes_amount+1):
    # Чтение файла с геном
    gene_path = f"{genomes_dir}/{species_dir}/Gene{file_num}.fasta"
    try:
        gene_record = next(SeqIO.parse(gene_path, "fasta"))
    except:
        raise ValueError(f"Gene reading erroe in {gene_path}")
    print(f"Gene {gene_path} is loaded")

    # Поиск и запись результатов
    find_and_write_fragment_positions(gene_record, fragment_records)
