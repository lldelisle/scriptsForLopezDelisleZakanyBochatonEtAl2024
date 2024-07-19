mkdir -p /scratch/ldelisle/HoxBstudy/baredSC
cd /scratch/ldelisle/HoxBstudy/baredSC

# We generate the table for parallel mcmc in 1d:
gitHubDirectory=/home/ldelisle/softwares/scriptsForLopezDelisleZakanyBochatonEtAl2024/
pathForTable="${gitHubDirectory}/scRNAseq/baredSC/table_baredSC_1d.txt"
echo -e "${gitHubDirectory}/scRNAseq/plots/meta_data_mut_NMPs.txt\tHoxb13\t3" >> $pathForTable

sbatch --array 1-1 --chdir $PWD/ ${gitHubDirectory}/scRNAseq/baredSC/sbatch_baredSC_1d.sh ${pathForTable} $PWD/baredSC_1d/

# Then we generate the table for parallel mcmc in 2d:
gitHubDirectory=/home/ldelisle/softwares/scriptsForLopezDelisleZakanyBochatonEtAl2024/
pathForTable="${gitHubDirectory}/scRNAseq/baredSC/table_baredSC_2d.txt"
for gene in "Hoxa9" "Hoxb9" "Hoxc9" "Hoxd9"; do
    echo -e "${gitHubDirectory}/scRNAseq/plots/meta_data_mut_NMPs.txt\t${gene}\tHoxb13\t3\t3" >> $pathForTable
done

sbatch --array 1-4 --chdir $PWD/ ${gitHubDirectory}/scRNAseq/baredSC/sbatch_baredSC_2d.sh ${pathForTable} $PWD/baredSC_2d/

tar zcvmf baredscToCheck.tar.gz */*convergence.png */*corner.png

tar zcvmf baredscZakanyBochaton20240326npz.tar.gz */*npz
mv baredscZakanyBochaton20240326npz.tar.gz /work/updub/
