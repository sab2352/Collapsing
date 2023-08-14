## Navigation
### Go back to <a href = "/Tutorial/Home.md">Wiki</a>
<hr>

# Selecting Models Options 
<h2> Here are the different parameters for each model and what they mean: </h2>

| ATAV Options | Notes |
| --- | ----------- |
| --collapsing-lite | This is a collapsing done using the Master Genotype File |
| --mann-whitney-test | We are performing a Mann-Whitney Test for collapsing |
| --genotype | This is your input. It is in the $genos variable. |
| --effect | This is how we pick the type of variant we are running the collapsing one (ie. Missense, PTV, Missense/PTV). Check out <a href = "#effects">this section</a> for what each categorization is. In general:<li> $MISSENSE_ONLY for missense only </li> <li> $FUNCTIONAL_EFFECTS for missense + PTV </li> <li> $FUNCTIONAL_EFFECTS_HQ is more stringent </li> <li> $LOF_EFFECTS for PTV </li> <li> $LOF_EFFECTS_HQ is more stringent </li> |
| --ensemble-missense | This picks missense variants that passed 2/3 of the Polyphen, Primate AI and Revel score criterias. |
| --polyphen | This is one out of three protein prediction scores for missense variants. Ensemble typically sets this to probably. A more lenient filter is probably,possibly,unknown. Can also be used without the --ensemble-missense option, which means that all 3 (PolyPhen, Primate-AI and Revel scores) have to be passed rather than just 2/3. |
| --min-primate-ai | This is one out of three protein prediction scores for missense variants. Ensemble typically sets this to 0.8. A more lenient filter is 0.5. Can also be used without the --ensemble-missense option, which means that all 3 (PolyPhen, Primate-AI and Revel scores) have to be passed rather than just 2/3. |
| --min-revel-score | This is one out of three protein prediction scores for missense variants. Ensemble typically sets this to 0.5. A more lenient filter is 0.35. Can also be used without the --ensemble-missense option, which means that all 3 (PolyPhen, Primate-AI and Revel scores) have to be passed rather than just 2/3. |
| --min-exac-vqslod-snv<br>--min--exac-vqslod-indel | The VQSLOD score is the log odds that a variant is a true variant as opposed to an artifact. Usually set to 5000. |
| --exclude-repeat-regions | Repeat regions are low-complexity regions. We want to exclude these because the allele frequencies regions for these regions are hard to calculate due to variations in repeat regioins. |
| --gnomad-genome-pop<br>--gnomad-exome-pop<br>--exac-pop | This is the population that the allele frequency is calculated from. 'global' means that the allele frequency is calculated from all the populations. Some other populations are: afr,amr,asj,eas,sas,fin,nfe |
| --gnomad-genome-maf<br>--gnomad-exome-maf<br>--exac-maf | This is the allele frequency of your model. |
| --loo-maf | Inside an ancestry-based cluster, the same variant cannot be present across multiple samples within the same ancesry, or it'll mean that the variant is common for that particular ancestry. In other words, this filter makes sure that the variant is actually rare. We want a filter that finds only 1 out of the ancestry to contain the variant. An example LOO (leave-one-out) value can be 0.00001.  |
| --max-qc-fail-sample | Usually set to 0 for ultra rare and 2 for rare.  |
| --min-pext-ratio | Pext ratio score from GnomAD is the proportion expressed across transcripts. It is usually set to 0.9. A 0.9 ratio means that the region is transcribed at least 90%. |
| &#x2011;&#x2011;max&#x2011;sub&#x2011;rvis&#x2011;domain&#x2011;score&#x2011;percentile | Usually set to 50. Tells you the degree of intolerance for the exon or protein domain in which a variant falls. |
| --exclude-false-loftee | Loss-Of-Function Transcript Effect Estimator (LOFTEE) determines the effects of loss-of-function variants. For PTV models.  |

<br>

<h2 id = "effects"> Effects </h2>

| Effects | Variants |
| --- | ----------- |
| CODING_SPLICE | HIGH:exon_loss_variant<br>HIGH:frameshift_variant<br>HIGH:rare_amino_acid_variant<br>HIGH:stop_gained<br>HIGH:start_lost<br>HIGH:stop_lost<br>HIGH:splice_acceptor_variant<br>HIGH:splice_donor_variant<br>HIGH:gene_fusion<br>HIGH:bidirectional_gene_fusion<br>MODERATE:3_prime_UTR_truncation+exon_loss_variant<br>MODERATE:5_prime_UTR_truncation+exon_loss_variant<br>MODERATE:coding_sequence_variant<br>MODERATE:disruptive_inframe_deletion<br>MODERATE:disruptive_inframe_insertion<br>MODERATE:conservative_inframe_deletion<br>MODERATE:conservative_inframe_insertion<br>MODERATE:missense_variant+splice_region_variant<br>MODERATE:missense_variant<br>MODERATE:splice_region_variant<br>LOW:5_prime_UTR_premature_start_codon_gain_variant<br>LOW:initiator_codon_variant<br>LOW:initiator_codon_variant+non_canonical_start_codon<br>LOW:splice_region_variant+synonymous_variant<br>LOW:splice_region_variant<br>LOW:start_retained<br>LOW:stop_retained_variant<br>LOW:synonymous_variant |
| FUNCTIONAL_EFFECTS | HIGH:exon_loss_variant<br>HIGH:frameshift_variant<br>HIGH:rare_amino_acid_variant<br>HIGH:stop_gained<br>HIGH:start_lost<br>HIGH:stop_lost,HIGH:splice_acceptor_variant<br>HIGH:splice_donor_variant<br>HIGH:gene_fusion<br>HIGH:bidirectional_gene_fusion<br>MODERATE:3_prime_UTR_truncation+exon_loss_variant<br>MODERATE:5_prime_UTR_truncation+exon_loss_variant<br>MODERATE:coding_sequence_variant<br>MODERATE:disruptive_inframe_deletion<br>MODERATE:disruptive_inframe_insertion<br>MODERATE:conservative_inframe_deletion<br>MODERATE:conservative_inframe_insertion<br>MODERATE:missense_variant+splice_region_variant<br>MODERATE:missense_variant<br>LOW:5_prime_UTR_premature_start_codon_gain_variant<br>LOW:initiator_codon_variant<br>LOW:initiator_codon_variant+non_canonical_start_codon |
| MOD_HIGH_EFFECTS | HIGH:exon_loss_variant<br>HIGH:frameshift_variant<br>HIGH:rare_amino_acid_variant<br>HIGH:stop_gained<br>HIGH:start_lost<br>HIGH:stop_lost<br>HIGH:splice_acceptor_variant<br>HIGH:splice_donor_variant<br>HIGH:gene_fusion<br>HIGH:bidirectional_gene_fusion<br>MODERATE:3_prime_UTR_truncation+exon_loss_variant<br>MODERATE:5_prime_UTR_truncation+exon_loss_variant<br>MODERATE:coding_sequence_variant<br>MODERATE:disruptive_inframe_deletion<br>MODERATE:disruptive_inframe_insertion<br>MODERATE:conservative_inframe_deletion<br>MODERATE:conservative_inframe_insertion<br>MODERATE:missense_variant+splice_region_variant<br>MODERATE:missense_variant<br>MODERATE:splice_region_variant |
| FUNCTIONAL_EFFECTS_HQ | HIGH:frameshift_variant<br>HIGH:stop_gained<br>HIGH:splice_acceptor_variant<br>HIGH:splice_donor_variant<br>MODERATE:missense_variant+splice_region_variant<br>MODERATE:missense_variant |
| LOF_EFFECTS | HIGH:exon_loss_variant<br>HIGH:frameshift_variant<br>HIGH:rare_amino_acid_variant<br>HIGH:stop_gained<br>HIGH:stop_lost<br>HIGH:start_lost<br>HIGH:gene_fusion<br>HIGH:bidirectional_gene_fusion<br>HIGH:splice_acceptor_variant<br>HIGH:splice_donor_variant |
| LOF_EFFECTS_HQ | HIGH:frameshift_variant<br>HIGH:stop_gained<br>HIGH:splice_acceptor_variant<br>HIGH:splice_donor_variant |
| SYN_EFFECTS | LOW:start_retained<br>LOW:stop_retained_variant<br>LOW:synonymous_variant |
| MISSENSE_ONLY | MODERATE:missense_variant+splice_region_variant<br>MODERATE:missense_variant |

<h2> Sample Models</h2>
<h3> Dominant Synonymous </h3>

<code>model_7: "$atav --collapsing-lite --mann-whitney-test --genotype $genos --effect $SYN_EFFECTS --min-exac-vqslod-snv 5000 --min-exac-vqslod-indel 5000  --gnomad-genome-pop global --gnomad-genome-maf 0 --gnomad-exome-pop global --gnomad-exome-maf 0 --exac-pop global --exac-maf 0 --loo-maf 0.0005 --min-pext-ratio 0.9 --max-qc-fail-sample 0 --exclude-repeat-region --sample $samples --out $outputFolder/domSynP9/"</code>


<h3> Ultra-rare PTV + Missense </h3>

<code>model_13: "$atav --collapsing-lite --mann-whitney-test --genotype $genos --effect $FUNCTIONAL_EFFECTS_HQ --polyphen probably --min-primate-ai 0.8 --min-revel-score 0.5 --ensemble-missense --min-exac-vqslod-snv 5000 --min-exac-vqslod-indel 5000  --gnomad-genome-pop global --gnomad-genome-maf 0 --gnomad-exome-pop global --gnomad-exome-maf 0 --exac-pop global --exac-maf 0 --loo-maf 0.0005 --max-qc-fail-sample 0 --exclude-repeat-region --exclude-false-loftee --min-pext-ratio 0.9 --sample $samples --out $outputFolder/URPTVMisP9/"</code>


<h3> Ultra-rare Missense </h3>

<code>model_14: "$atav --collapsing-lite --mann-whitney-test --genotype $genos --effect $MISSENSE_ONLY --polyphen probably --min-primate-ai 0.8 --min-revel-score 0.5 --ensemble-missense --min-exac-vqslod-snv 5000 --min-exac-vqslod-indel 5000  --gnomad-genome-pop global --gnomad-genome-maf 0 --gnomad-exome-pop global --gnomad-exome-maf 0 --exac-pop global --exac-maf 0 --loo-maf 0.0005 --max-qc-fail-sample 0 --min-pext-ratio 0.9 --exclude-repeat-region --sample $samples --out $outputFolder/URMisP9/"</code>


<h3> Dominant Ultra-rare PTV </h3>

<code>model_10: "$atav --collapsing-lite --mann-whitney-test --genotype $genos --effect $LOF_EFFECTS_HQ --min-exac-vqslod-snv 5000 --min-exac-vqslod-indel 5000  --gnomad-genome-pop global --gnomad-genome-maf 0 --gnomad-exome-pop global --gnomad-exome-maf 0 --exac-pop global --exac-maf 0 --loo-maf 0.0005 --max-qc-fail-sample 0  --exclude-false-loftee --min-pext-ratio 0.9 --exclude-repeat-region --sample $samples --out $outputFolder/URPTVP9/"</code>

