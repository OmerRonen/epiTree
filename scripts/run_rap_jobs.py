import os.path
import subprocess

import numpy as np
import pandas as pd


def save_files(phenotype):
    """
    Get the individuals for a given phenotype
    """
    pth = os.path.join(f"data/{phenotype}")
    if not os.path.exists(pth):
        os.makedirs(pth)
    data_pheno = pd.read_csv("disease_phenotypes_full_PC40.csv")
    data_pheno = data_pheno.loc[:, ["participant.eid", phenotype]]
    # take all positive cases and sample the same number of controls at random
    data_pheno_case = data_pheno[data_pheno[phenotype] == 1]
    data_pheno_control = data_pheno[data_pheno[phenotype] == 0].sample(n=data_pheno_case.shape[0])
    data_pheno = np.concatenate([data_pheno_control, data_pheno_case], axis=0)

    # read fam file as tsv
    fam_data = np.genfromtxt("epitree/pick.fam", delimiter='\t')
    idx = [id in data_pheno[:, 0] for id in fam_data[:, 0]]
    fam_data = fam_data[idx, :].astype(int)

    # save the fam data as tsv
    np.savetxt(os.path.join(pth, "ind.fam"), fam_data, fmt='%s', delimiter='\t')

    # save the data_pheno as csv
    data_pheno = pd.DataFrame(data_pheno).astype(int)
    data_pheno.to_csv(os.path.join(pth, f"pheno.csv"), index=False, header=False)
    individuals = data_pheno.iloc[:, 0].values
    # save into a mapping file with to columns of individual id and individual id
    df = pd.DataFrame({"1":individuals, "2":individuals})
    df.to_csv(os.path.join(pth, f"mapping.csv"), index=False, header=False)

    # # save as txt file each individual as a row
    # np.savetxt(os.path.join(pth, f"ind.txt"), individuals, fmt='%s')


def run_job(chr, tissue, phenotype):

    # save_files(phenotype)
    # cmd_upload = f"dx upload -r data/{phenotype}"
    # # run and wait for the upload to finish
    # subprocess.run("conda activate /scratch/users/omer_ronen/mutemb", shell=True)
    # subprocess.run(cmd_upload, shell=True)
    phenotype = phenotype.replace(" ", "_")
    pth_data = f"{phenotype}"
    # tissue = "Brain_Cortex"

    geno_path = f"project-GJkXqkQJ5Xvyq6KQ708Q3qqG:/Bulk/Imputation/UKB imputation from genotype"
    sample_file = f"ukb22828_c{chr}_b0_v3.sample"
    bgen_file = f"ukb22828_c{chr}_b0_v3.bgen"
    snplist_file = f"snpList_ctimp_{tissue}.snplist"
    individuals_file = "ind.fam"
    out_prefix = f"{tissue}_{phenotype}_chr{chr}".replace(" ", "_")
    files_str = f'-iin="{os.path.join(geno_path,sample_file)}" -iin="{os.path.join(geno_path,bgen_file)}" -iin="{snplist_file}" -iin="{os.path.join(pth_data, "ind.fam")}"'
    plink_cmd = f'"plink2 --bgen {bgen_file} ref-first --sample {sample_file} --extract {snplist_file} --keep {individuals_file} --make-bed --out {out_prefix}"'
    cmd = f'/scratch/users/omer_ronen/mutemb/bin/dx run --instance-type mem1_ssd1_v2_x36 --yes swiss-army-knife {files_str} -icmd={plink_cmd}'
    # execute the command
    subprocess.run(cmd,  shell=True)
    print(cmd)


if __name__ == '__main__':
    # save_files("Systemic lupus erythematosus")
    # save_files("Rheumatoid arthritis")
    pheno = "Rheumatoid arthritis"
    expression = "Whole_Blood"
    # for chr in range(1, 22):
    run_job(22,expression ,"Rheumatoid arthritis")
    run_job(22,expression ,"Systemic lupus erythematosus")

    # run_job(3, "Multiple sclerosis")