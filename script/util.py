import pysam


def load_aln(bam_file, aln_file):
    with pysam.AlignmentFile(bam_file, 'rb') as fin, open(aln_file, 'wt') as fout:
        fout.write('q_name\tq_size\tq_start\tq_end\tt_name\tt_size\tt_start\t'
                   't_end\tstrand\tlocal_identity\tglobal_identity\tcigar\tq_seq\n')

        for record in fin.fetch():
            if record.is_unmapped:
                continue
            strand = '-' if record.is_reverse else '+'
            q_name = record.query_name
            q_size = record.query_length
            cigartuples = record.cigartuples
            q_start = 1
            cigartuple0 = cigartuples[0]
            if strand == '-':
                cigartuple0 = cigartuples[-1]
            op0 = cigartuple0[0]
            op_len0 = cigartuple0[1]
            if op0 == 4:
                q_start = op_len0 + 1
            m_num = 0
            i_num = 0
            d_num = 0
            s_num = 0
            x_num = 0
            for cigartuple in cigartuples:
                op = cigartuple[0]
                op_len = cigartuple[1]

                if op == 0:
                    m_num += op_len
                elif op == 1:
                    i_num += op_len
                elif op == 2:
                    d_num += op_len
                elif op == 4:
                    s_num += op_len
                else:
                    x_num += op_len
            if x_num > 0:
                continue
            t_name = record.reference_name
            t_size = record.header.get_reference_length(record.reference_name)
            q_align_len = m_num + i_num
            q_end = q_start + q_align_len - 1
            t_start = record.reference_start
            t_align_len = m_num + d_num
            t_end = t_start + t_align_len - 1
            ed = record.get_tag("NM")
            local_identity = (q_align_len - ed) / q_align_len
            diff = s_num + i_num + d_num + ed
            global_identity = (q_size - diff) / q_size
            fout.write(f'{q_name}\t{q_size}\t{q_start}\t{q_end}\t'
                       f'{t_name}\t{t_size}\t{t_start}\t{t_end}\t'
                       f'{strand}\t{local_identity:.2f}\t{global_identity:.2f}\t'
                       f'{record.cigarstring}\t{record.query_sequence}\n')


def filter_fq(read_names_file, fq_file, anchored_fq_file):
    read_names = set()
    with open(read_names_file, 'rt') as fin:
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            read_names.add(line)

    with pysam.FastqFile(fq_file) as fin, open(anchored_fq_file, 'wt') as fout:
        for entry in fin:
            name = entry.name
            start = name.find(' ')
            if start != -1:
                name = name[0:start]
            if name in read_names:
                fout.write('@' + name + '\n')
                fout.write(entry.sequence + '\n')
                fout.write("+" + '\n')
                fout.write(entry.quality + '\n')



def filter_sam_header(chr_name, in_file, fout):
    fin = open(in_file, 'rt')
    for line in fin:
        line = line.strip()
        if len(line) == 0:
            continue

        ss = line.split('\t')
        if line.startswith('@SQ'):
            if ss[1] == 'SN:' + chr_name:
                fout.write(line + '\n')
    fin.close()


def filter_sam(read_name, in_file, fout):
    fin = open(in_file, 'rt')
    for line in fin:
        line = line.strip()
        if len(line) == 0:
            continue

        ss = line.split('\t')
        if read_name == ss[0]:
            fout.write(line + '\n')

    fin.close()
