# -*- coding: UTF-8 -*-
import os
import sys
import shutil
import ConfigParser
from helpers import *




def main():
# ------------------------------------------------------------------------- #
#       1. CONFIGURATION & OUTPUT SET UP                                    #
# ------------------------------------------------------------------------- #
    
    config = ConfigParser.ConfigParser()
    config.readfp(open(r'config.txt'))

    prokka = config.get('configuration', 'prokka')
    blastn = config.get('configuration', 'blastn')
    blastp = config.get('configuration', 'blastp')
    makedb = config.get('configuration', 'makeblastdb')
    ccr = config.get('configuration', 'ccr')
    orfX = config.get('configuration', 'orfX')
    pbp2a = config.get('configuration', 'pbp2a')
    attr_db = config.get('configuration', 'attr_db')
    input_dir = config.get('configuration', 'input_dir')
    output_data = config.get('configuration', 'output')
    threshold = config.get('configuration', 'threshold')
    core_db = config.get('configuration', 'core_db')
    base_network = config.get('configuration', 'base_network')
    chunk = config.get('configuration', 'chunk')
    mash = config.get('configuration', 'mash')
    core_dir = config.get('configuration', 'core_dir')
    core_proteins_sccmec = config.get('configuration', 
        'core_proteins_sccmec')
    features = config.get('configuration', 'features')
    labels = config.get('configuration', 'labels')

    basedir = os.path.dirname(os.path.abspath(sys.argv[0]))
    
    mash_distances = os.path.join(basedir, features)
    sccmectypes = os.path.join(basedir, labels)

    core_dir_path = os.path.join(basedir, core_dir)
    core_dict_file = os.path.join(basedir, core_proteins_sccmec)
    core_proteins_sccmec = {}
    with open(core_dict_file) as fileobj:
        core_proteins_sccmec = dict(line.strip().split(',', 1) 
            for line in fileobj)

    mash_path = os.path.join(basedir, mash)
    makeblastdb_exe = os.path.join(basedir, makedb)
    prokka_exe = blastp_exe = os.path.join(basedir, prokka)
    blastn_exe = os.path.join(basedir, blastn)
    blastp_exe = os.path.join(basedir, blastp)
    input_contigs = os.path.join(basedir, input_dir)
    orfx_base = os.path.join(basedir, orfX)
    ccr_base = os.path.join(basedir, ccr)
    pbp2a_base = os.path.join(basedir, pbp2a)
    global_output_dir = create_dir(basedir, output_data)
    attr_database_path = os.path.join(basedir, attr_db)
    core_db_path = os.path.join(basedir, core_db)
    base_network_path = os.path.join(basedir, base_network)
    chunk_path = os.path.join(basedir, chunk)


# ------------------------------------------------------------------------- #
#       2. LOOP THROUGH CONTIGS FILES IN INPUT_DIR                          #
# ------------------------------------------------------------------------- #

    contigs_to_analize = [f for f in os.listdir(input_contigs) if 
                            f.endswith('.fasta')]

    for contig in contigs_to_analize:
# ------------------------------------------------------------------------- #
#       a. Create output folder for contig output_<name>                    #
        print(contig)
        contig_abspath = os.path.join(input_contigs, contig)
        filename = contig.split('.')[0]
        folder = 'output_{0}'.format(filename)
        output_folder = create_dir(global_output_dir, folder)

        os.chdir(output_folder)


# ------------------------------------------------------------------------- #
#       b. Run prokka in contig file                                        # 
#       keep fna and ffn files in dicts to be used to determine             #   
#       in which contigs are codified orfx, mec and ccrs.                   #

        output_prokka = 'output_prokka_{0}'.format(filename)
        execute_prokka(prokka_exe, output_prokka, contig_abspath, filename)
        ffn, gff, fna, faa = prokka_files(output_prokka)
        ffn, gff, fna, faa = list(os.path.join(output_folder, x) 
            for x in [ffn, gff, fna, faa])
        fna_dict = fasta2dict(fna)
        ffn_dict = fasta2dict(ffn)

# ------------------------------------------------------------------------- #
#       c. BLAST databases configuration                                    # 
#       To perform blastn searches for orfx gene                            #
#       and blastp searches for mec and ccrs (different filter parameters)  # 

        nucl_db_dir = create_dir(output_folder, 'nucl_db_dir')
        nucl_db = makeblastdb(makeblastdb_exe, nucl_db_dir, ffn, 
            'nucl', 'nucl_db')
        nucl_db_path = os.path.join(nucl_db_dir, nucl_db)

        prot_db_dir = create_dir(output_folder, 'prot_db_dir')
        prot_db = makeblastdb(makeblastdb_exe, prot_db_dir, faa, 
            'prot', 'prot_db')
        prot_db_path = os.path.join(prot_db_dir, prot_db)

# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
#       d. Check if is MRSA by running BLAST searches of core elements      #
#                                                                           #
#       FIRST data validation                                               # 
#       -> mrsa if and only if ccr and mec and orfx == True                 #
#       changes can be made to ignore abscence of orfx.                     # 
#                                                                           #
#       Continue with next contig to analize if at least one core element   #
#       is not present.                                                     #

        orfx_nucl_hit, err = simpleBlast(blastn_exe, nucl_db_path, 
            orfx_base, 'orfX')
        orfx_nucl_hit = simpleBlastParser(orfx_nucl_hit)
        mec_hit, err = simpleBlast(blastp_exe, prot_db_path, 
            pbp2a_base, 'PBP2a')
        mec_hit = blast_mec_parser(mec_hit)
        ccr_hit, err = simpleBlast(blastp_exe, prot_db_path, 
            ccr_base, 'ccr')
        ccr_hit = blast_ccr_parser(ccr_hit)

        print('orfX Hit: ', orfx_nucl_hit)
        print('PBP2a Hit: ', mec_hit)
        print('ccr Hit: ', ccr_hit)

        print(nucl_db_dir)
        print(prot_db_dir)
        prokkk = os.path.join(output_folder, output_prokka)
        print(prokkk)
        core_elements = [orfx_nucl_hit, mec_hit, ccr_hit]

        if not all(core_elements):
            print('It\'s not MRSA')

            shutil.rmtree(nucl_db_dir, ignore_errors=False, onerror=None)
            shutil.rmtree(prot_db_dir, ignore_errors=False, onerror=None)
            shutil.rmtree(prokkk, ignore_errors=False, onerror=None)
            os.chdir(output_folder)
            info_file = 'MSSA_{0}.txt'.format(filename)
            with open(info_file, "w") as f:
                f.write('Gene\tLocus\n')
                f.write('orfX\t{0}\n'.format(str(core_elements[0])))
                f.write('PBP2a\t{0}\n'.format(str(core_elements[1])))
                f.write('ccr\t{0}\n'.format(str(core_elements[2])))
            if core_elements[0] is not None:
                actual_orfx = get_sequence(ffn_dict, orfx_nucl_hit)
                attB_sequence = actual_orfx[len(actual_orfx)-60:]
                attB_datafile = 'attB_{0}.fasta'.format(filename)
                with open(attB_datafile, "w") as f:
                    f.write('>attB_{0}\n'.format(filename))
                    f.write(attB_sequence+'\n')

            continue

# ------------------------------------------------------------------------- #
#              PRIMERA PARTE: EXTRAER SCCMEC DESDE LOS CONTIGS              #
# ------------------------------------------------------------------------- #

# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
#       e. Check if CORE ELEMENTS are in the same contig                    #
#       after check if it is MRSA, find out in which contigs orfx, mec and  #
#       ccr are present. Here we are looking for the contig header of each  #
#       gene to be compared in check_region_id()                            #

        print('It\'s MRSA'+'\n')
        os.chdir(output_folder)

        actual_orfx = get_sequence(ffn_dict, orfx_nucl_hit)
        actual_mec = get_sequence(ffn_dict, mec_hit)
        actual_ccr = get_sequence(ffn_dict, ccr_hit)

        seq_region_id_orfx = get_region(gff, orfx_nucl_hit)
        seq_region_id_pbp2a = get_region(gff, mec_hit)
        seq_region_id_ccr = get_region(gff, ccr_hit)

        region_ids = [seq_region_id_orfx, 
            seq_region_id_pbp2a, 
            seq_region_id_ccr]

        print('orfX in: ', seq_region_id_orfx, 
            len(get_sequence(fna_dict, seq_region_id_orfx)))
        print('PBP2a in: ', seq_region_id_pbp2a, 
            len(get_sequence(fna_dict, seq_region_id_pbp2a)))
        print('ccr in: ', seq_region_id_ccr, 
            len(get_sequence(fna_dict, seq_region_id_ccr)))


        cassette_filename = 'sccmec_{0}.fasta'.format(filename)

        if check_region_id(region_ids):
# ------------------------------------------------------------------------- #
#       f. If all core elements are the same fna sequence                   #
#                                                                           #
#       checkSense() parses the fna sequences in which orfx is present      # 
#       to keep attL to the left side (our point of view.)                  #

            fna_sequence = get_sequence(fna_dict, seq_region_id_orfx)
            fna_sequence = checkSense(actual_orfx, fna_sequence)
            print('length contig with orfX: ', len(fna_sequence))

#       We are only interested in the sequences from orfx to the end        #
#       template_dna will contain that given subsequences                   #
#       and it will be used to search the right end of the SCCmec           #
            texto_orfx, inicio_orfx, final_orfx = get_orfX_pos(fna_sequence, 
                actual_orfx, 'orfX')
            print('orfX info: ', texto_orfx, inicio_orfx, final_orfx)
            template_dna = fna_sequence[inicio_orfx:]

# ------------------------------------------------------------------------- #
#       G. Find attachment site sequences right (attR)                      #
#       As we know, attR is present at the right end of the SCCmec          #
#       this ~20pb sequences contains a core pattern of 8 nucleotides, this #
#       pattern is being search using re and an extended hit (58nt) is used #
#       as a query againts our attR database using blastn to probe that     # 
#       is a real attR.                                                     #
#       get_attchment (lol a typo) function computes the search             #
#       and return the location of the hit.                                 #

            params = get_attchment(filename, 
                                output_folder, 
                                template_dna, 
                                attr_database_path, 
                                blastn_exe)

            if params:
                attr_s, attr_e, hit = params
                orfx_text, orfx_s, orfx_e = sequence_position(template_dna, 
                    actual_orfx, "orfX")
                att_actual_orfx = actual_orfx[len(actual_orfx)-60:]
                attL_datafile = 'attL_{0}.fasta'.format(filename)
                with open(attL_datafile, 'w') as f:
                    f.write('>attL_{0}_{1}_{2}\n'.format(filename, 
                        orfx_s, orfx_e))
                    f.write(att_actual_orfx+'\n')

#       cassette variable will contain the sequences from the               # 
#       beginning of orfX to the attR-end.                                  # 
                cassette = template_dna[orfx_s:attr_e]

#################################-BUG-PRONE-#################################
                # Este es un Bug que encontre provando sequencias al azar
                #y es que se necesita comprobar que la secuencia en la
                # variable "cassette" posee los tres elementos
                # creo que lo obvie porque tenemos cerca de 200 secuencias
                # attR, pero si aparece un cassette muy divergente vamos
                # a tener problemas
                try: 
                    mec_text, mec_s, mec_e = sequence_position(cassette, 
                        actual_mec, "PBP2a")
                except TypeError:
                    info_file = 'MRSA_{0}_wo_mecA.txt'.format(filename)
                    with open(info_file, 'w') as f:
                        f.write('mecA out of selected range\n')
                    shutil.rmtree(nucl_db_dir, ignore_errors=False, onerror=None)
                    shutil.rmtree(prot_db_dir, ignore_errors=False, onerror=None)
                    shutil.rmtree(prokkk, ignore_errors=False, onerror=None)
                    continue

                try:
                    ccr_text, ccr_s, ccr_e = sequence_position(cassette, 
                        actual_ccr, "ccr")
                except TypeError:
                    info_file = 'MRSA_{0}_wo_ccr.txt'.format(filename)
                    with open(info_file, 'w') as f:
                        f.write('ccr out of selected range\n')
                    shutil.rmtree(nucl_db_dir, ignore_errors=False, onerror=None)
                    shutil.rmtree(prot_db_dir, ignore_errors=False, onerror=None)
                    shutil.rmtree(prokkk, ignore_errors=False, onerror=None)
                    continue
#############################################################################

# ------------------------------------------------------------------------- #
#       H. Create a file with SCCmec sequences from orfX to attR            #
#       The extraction of SCCmec is qualify as A given that all elements    #
#       (orfX, mec and crr) were in the same contig and there wasn't        #
#       problems to find the attR.                                          #   
# ------------------------------------------------------------------------- #
                 

                with open(cassette_filename, 'w') as f:
                    f.write('>sccmec_{0}_l{1}\n'.format(filename, 
                        str(len(cassette))))
                    for i in range(0, len(cassette), 60):
                        f.write(cassette[i:i+60]+'\n')

                info_file = 'MRSA_{0}.txt'.format(filename)

                with open(info_file, 'w') as f:
                    f.write('fna length\t{0}\n'.format(
                        str(len(template_dna))))
                    f.write('SCCmec length\t{0}\n'.format(
                        str(len(cassette))))                
                    f.write('Quality\tA\n')


                cassette_path = os.path.join(output_folder, 
                    cassette_filename)
        
# ------------------------------------------------------------------------- #
#       i. If attR is not found in the same contig where core elemets were  #
#       found, the script end the search and continue with the next file.   #
#       (very unlikely scenario)                                            #
# ------------------------------------------------------------------------- #
            else:
                shutil.rmtree(nucl_db_dir, ignore_errors=False, onerror=None)
                shutil.rmtree(prot_db_dir, ignore_errors=False, onerror=None)
                shutil.rmtree(prokkk, ignore_errors=False, onerror=None)
                
                info_file = 'attr_not_found_{0}.txt'.format(filename)
                with open(info_file, 'w') as f:
                    f.write('attr not found\n')
                continue


# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
#       J. If each fna contain a core elemets the script continue with the  #
#       next file to analize. In the case of two differents fna read K.     #
#       (note that i use fna and contig in a interchangeably way and i hope #
#       it doesn't confuse whoever review this code)                        #
        else:

            print('core elements contained in more than one fna')

            fnas = sort_fna(seq_region_id_orfx, 
                seq_region_id_pbp2a, 
                seq_region_id_ccr)
            info_file = 'separated_core_elements_{0}.txt'.format(filename)
            if len(fnas) == 3:
                shutil.rmtree(nucl_db_dir, ignore_errors=False, onerror=None)
                shutil.rmtree(prot_db_dir, ignore_errors=False, onerror=None)
                shutil.rmtree(prokkk, ignore_errors=False, onerror=None)
                
                with open(info_file, 'w') as f:
                    f.write('each core elements in differents fna\n')
                continue

# ------------------------------------------------------------------------- #
#       K. To emulate what was made in G, we first need to find which fna   #
#       contains: orfX+mec, orfX+ccr or just orfX. This is being made by    #
#       sort_fna() function (line ~277). By doing this we will be able to   #
#       find the left end (orfX presents here) and we can look for the attR #
#       in the remaining fna sequence. The tricky part is that we don't     #
#       know if attR will appear in the same sense of the fna, thus we are  #
#       going to test both senses and we will only keep the sense in which  #
#       the attR was found.                                                 #            
            else:
                os.chdir(output_folder)

                contig_left, contig_right = fnas
                sequence_contig_left = get_sequence(fna_dict, contig_left)
                sequence_contig_right = get_sequence(fna_dict, contig_right)

                sequence_contig_left = checkSense(actual_orfx, 
                    sequence_contig_left)
                orfx_text, orfx_s, orfx_e = sequence_position(
                    sequence_contig_left, actual_orfx, "orfX")
                
                att_actual_orfx = actual_orfx[len(actual_orfx)-60:]
                attL_datafile = 'attL_{0}.fasta'.format(filename)
                with open(attL_datafile, 'w') as f:
                    f.write('>attL_{0}_{1}_{2}\n'.format(filename, 
                        orfx_s, orfx_e))
                    f.write(att_actual_orfx+'\n')

                cassette_left_end = sequence_contig_left[orfx_s:]

# ------------------------------------------------------------------------- #
#       L. Search for attR using the same sense                             #
                params = get_attchment(filename, 
                                    output_folder, 
                                    sequence_contig_right, 
                                    attr_database_path, 
                                    blastn_exe)
                if params:
                    print('get attR 1')
                    attr_s, attr_e, hit = params
                    print(hit)
                    cassette_right_end = sequence_contig_right[:attr_e]

                    cassette = ''.join([cassette_left_end, 
                        cassette_right_end])
#################################-BUG-PRONE-#################################
                    try: 
                        mec_text, mec_s, mec_e = sequence_position(cassette, 
                            actual_mec, "PBP2a")
                    except TypeError:
                        info_file = 'MRSA_{0}_wo_mecA.txt'.format(filename)
                        with open(info_file, 'w') as f:
                            f.write('mecA out of selected range\n')
                        shutil.rmtree(nucl_db_dir, ignore_errors=False, onerror=None)
                        shutil.rmtree(prot_db_dir, ignore_errors=False, onerror=None)
                        shutil.rmtree(prokkk, ignore_errors=False, onerror=None)
                        continue

                    try:
                        ccr_text, ccr_s, ccr_e = sequence_position(cassette, 
                            actual_ccr, "ccr")
                    except TypeError:
                        info_file = 'MRSA_{0}_wo_ccr.txt'.format(filename)
                        with open(info_file, 'w') as f:
                            f.write('ccr out of selected range\n')
                        shutil.rmtree(nucl_db_dir, ignore_errors=False, onerror=None)
                        shutil.rmtree(prot_db_dir, ignore_errors=False, onerror=None)
                        shutil.rmtree(prokkk, ignore_errors=False, onerror=None)
                        continue


                    with open(cassette_filename, 'w') as f:
                        f.write('>sccmec_{0}_l{1}\n'.format(filename, 
                            str(len(cassette))))
                        for i in range(0, len(cassette), 60):
                            f.write(cassette[i:i+60]+'\n')

                    info_file = 'MRSA_{0}.txt'.format(filename)

                    with open(info_file, 'w') as f:
                        f.write('fna-L length\t{0}\n'.format(
                            str(len(sequence_contig_left))))
                        f.write('fna-R length\t{0}\n'.format(
                            str(len(sequence_contig_right))))
                        f.write('SCCmec length\t{0}\n'.format(
                            str(len(cassette))))                
                        f.write('Quality\tB\n')

                    cassette_path = os.path.join(output_folder, 
                        cassette_filename)
                    
# ------------------------------------------------------------------------- #
#       M. Search for attR using the reverse complement sequences of fna    #
                else:
                    print('get attR 2')
                    sequence_contig_right = reverse_complement(
                        sequence_contig_right)
                    params = get_attchment(filename, 
                                        output_folder, 
                                        sequence_contig_right, 
                                        attr_database_path, 
                                        blastn_exe)
                    if params:
                        attr_s, attr_e, hit = params
                        print(hit)
                        cassette_right_end = sequence_contig_right[:attr_e]

                        cassette = ''.join([cassette_left_end, 
                            cassette_right_end])

#################################-BUG-PRONE-#################################

                        try: 
                            mec_text, mec_s, mec_e = sequence_position(cassette, 
                                actual_mec, "PBP2a")
                        except TypeError:
                            info_file = 'MRSA_{0}_wo_mecA.txt'.format(filename)
                            with open(info_file, 'w') as f:
                                f.write('mecA out of selected range\n')
                            shutil.rmtree(nucl_db_dir, ignore_errors=False, onerror=None)
                            shutil.rmtree(prot_db_dir, ignore_errors=False, onerror=None)
                            shutil.rmtree(prokkk, ignore_errors=False, onerror=None)
                            continue

                        try:
                            ccr_text, ccr_s, ccr_e = sequence_position(cassette, 
                                actual_ccr, "ccr")
                        except TypeError:
                            info_file = 'MRSA_{0}_wo_ccr.txt'.format(filename)
                            with open(info_file, 'w') as f:
                                f.write('ccr out of selected range\n')
                            shutil.rmtree(nucl_db_dir, ignore_errors=False, onerror=None)
                            shutil.rmtree(prot_db_dir, ignore_errors=False, onerror=None)
                            shutil.rmtree(prokkk, ignore_errors=False, onerror=None)
                            continue


                        with open(cassette_filename, 'w') as f:
                            f.write('>sccmec_{0}_l{1}\n'.format(filename, 
                                str(len(cassette))))
                            for i in range(0, len(cassette), 60):
                                f.write(cassette[i:i+60]+'\n')

                        info_file = 'MRSA_{0}.txt'.format(filename)
                        with open(info_file, 'w') as f:
                            f.write('fna-L length\t{0}\n'.format(
                                str(len(sequence_contig_left))))
                            f.write('fna-R length\t{0}\n'.format(
                                str(len(sequence_contig_right))))
                            f.write('SCCmec length\t{0}\n'.format(
                                str(len(cassette))))                
                            f.write('Quality\tC\n')

                        cassette_path = os.path.join(output_folder, 
                            cassette_filename)

# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
#       O. If attR is not found in either of both sense (5'-3' o 3'-5')     #
#       the script continue with the next file to analize                   #
                    else:
                        shutil.rmtree(nucl_db_dir, ignore_errors=False, onerror=None)
                        shutil.rmtree(prot_db_dir, ignore_errors=False, onerror=None)
                        shutil.rmtree(prokkk, ignore_errors=False, onerror=None)
                        info_file = 'attr_not_found_{0}.txt'.format(
                            filename)
                        with open(info_file, 'w') as f:
                            f.write('attr not found\n')
                        continue


# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
#              SEGUNDA PARTE: ANOTACION,  CLASIFICACION Y                   #
#                             VISUALIZACION DE SCCMEC                       #      
# ------------------------------------------------------------------------- #

        print(cassette_path)

        sccmec_id, faa_file_sccmec, annotation_file = current_annotation(
            cassette_path, 
            blastp_exe, 
            prokka_exe, 
            core_db_path, 
            output_folder)

        print(base_network_path)


        close_cassette, bool_type, predicted = network_classification(
            threshold, 
            base_network_path, 
            mash_path, 
            chunk_path, 
            cassette_path, 
            sccmec_id,
            mash_distances,
            sccmectypes,)


        length_cassette = len(cassette)
        print(close_cassette)
        core_annotation = core_proteins_sccmec.get(close_cassette)
        print(core_annotation)
        core_proteins_file, tipo = core_annotation.split(',')

        print(tipo, predicted)
        sccmec_type_file = 'sccmec_cla_type_{0}.txt'.format(sccmec_id) 

        with open(sccmec_type_file, 'w') as f:
            if bool_type:
                f.write('new-type\n')
                f.write('related to group {0}, specific type {1}\n'.format(tipo, predicted))
            else:
                f.write('add-type\n')
                f.write('related to group {0}, specific type {1}\n'.format(tipo, predicted))

        print(core_proteins_file)

        selected_core_proteins = os.path.join(core_dir_path, 
            core_proteins_file)


        vis_sccmec(faa_file_sccmec, annotation_file, 
            length_cassette, selected_core_proteins, blastp_exe)

        shutil.rmtree(nucl_db_dir, ignore_errors=False, onerror=None)
        shutil.rmtree(prot_db_dir, ignore_errors=False, onerror=None)
        shutil.rmtree(prokkk, ignore_errors=False, onerror=None)

        
    sys.exit('done')


if __name__ == '__main__':
    main()

