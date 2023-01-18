

```bash
SITE='2'
YEARS=('2018' '2019' '2020')
CHROMS=('1' '2' '3' '4' '5')
for year in ${YEARS[@]}; do
for chrom in ${CHROMS[@]}; do
    echo "sbatch scripts/averageAF/wrapper_averageAF.sh $SITE $year $chrom"
    sbatch scripts/averageAF/wrapper_averageAF.sh $SITE $year $chrom
done
done

```