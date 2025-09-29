**Investigating Non-Invasive Biomarkers and Risk Factors Associated with
Liver Cirrhosis**

**1. Introduction**

Cirrhosis refers to the scarring of liver tissue which occurs as a
result of long-term damage sustained by various chronic liver diseases
and conditions [^1,2^](https://www.zotero.org/google-docs/?bEWnj8).
Predominantly, these diseases include viral hepatitis (specifically HBV
& and HCV), alcoholic liver disease (ALD), and non-alcoholic fatty liver
disease (NAFLD) [^1,3,4^](https://www.zotero.org/google-docs/?oQhsuD).
Disease progression varies based on the underlying cause, the presence
or absence of treatment, and ongoing liver injury. In advanced stages,
liver cirrhosis can affect and compromise the proper functioning of the
liver due to long-term damage, with these effects potentially becoming
irreversible and life-threatening
[^4^](https://www.zotero.org/google-docs/?0U2qzg).

Although ranking as the 11th overall leading cause of death and 3rd
leading cause of death in people aged 45-64 years in 2021, liver
cirrhosis lacks an enhanced cure
[^5^](https://www.zotero.org/google-docs/?hQKV1h). Current treatments
are focused on symptomatic relief and slowing down disease progression
for individuals with mild to moderate severity. In severe cases, where
the liver is significantly damaged, liver transplant remains the primary
treatment option [^1^](https://www.zotero.org/google-docs/?ow2G7c).

Recognizing the need for innovative therapeutic approaches, there is an
urgent need to explore non-invasive biomarkers and identify risk factors
associated with liver cirrhosis to advance our comprehension of the
disease[^2^](https://www.zotero.org/google-docs/?Dc9l6F). In this
exploration lies the potential to unravel the intricate relationship
between liver health and the gut microbiome, a connection that has
recently been implicated in numerous
diseases[^1,3,6^](https://www.zotero.org/google-docs/?gwELIz).

Studies involving patients with different levels and types of liver
disease were conducted to observe the changes in the microbiota over
time and during disease progression. In 2011, Yan et. al presented their
findings in observing gut dysbiosis in mice with ALD using 16S rRNA
sequencing, highlighting intestinal bacterial overgrowth in diseased
mice [^7^](https://www.zotero.org/google-docs/?kmtiIm). In 2013, Bajaj
et. al investigated these changes using pyrosequencing techniques and
ribosomal taxa analysis and reported significant differences in the gut
microbiota of patients with cirrhosis compared to those without
[^8,9^](https://www.zotero.org/google-docs/?K0ntfI). Subsequently,
Loomba et. al (2017) set out to establish gut microbiota-derived
signatures to predict the presence of advanced fibrosis and
characterized the association between the severity of fibrosis in NAFLD
patients, and gut dysbiosis
[^10^](https://www.zotero.org/google-docs/?RrElqq). However, most
studies that exist in the existing literature rely on 16S rRNA
sequencing, offering limited insights compared to the holistic view
provided by metagenomic shotgun sequencing
[^11^](https://www.zotero.org/google-docs/?ZNHbkg).

To comprehensively investigate non-invasive biomarkers and identify risk
factors associated with liver disease, a thorough analysis of the gut
microbiome is imperative. This necessitates using extensive sequencing
techniques, such as whole metagenome shotgun sequencing, to capture the
complete microbial landscape and enable a more nuanced understanding of
the relationship between liver cirrhosis and the gut microbiome.

This report aims to explore non-invasive biomarkers and potential
factors associated with liver cirrhosis using the comprehensive data
obtained through whole metagenome shotgun sequencing as collected by Qin
et al (2014) [^2^](https://www.zotero.org/google-docs/?vHkPH7).

**2. Methods**

**2.1 Sample Data**

A total of 97 subjects (NC = 43, ND = 54) participated in this study.
The diagnosis of liver cirrhosis in affected subjects adhered to
international guidelines, including the requirement of a biopsy.
Borderline or inconclusive cases were excluded from the diseased
population[^2^](https://www.zotero.org/google-docs/?sf3ap7).
Additionally, individuals with hypertension, diabetes, metabolic
syndrome, inflammatory bowel disease (IBD), NAFLD, coeliac disease, and
cancer were excluded from the control group. Figure 1 presents
additional characteristics of the participants involved. Stool samples
from each subject were collected and a DNA library was constructed
according to Illumina instructions to perform whole metagenomic shotgun
sequencing. The same workflows from Illumina were then used to perform
whole metagenome shotgun sequencing, and reads that mapped to the human
genome, along with their mated/paired reads, were eliminated from each
sample using BWA with '-n 0.2'.

> ![](media/image1.png){width="3.9218755468066493in"
> height="2.9190660542432196in"}

**Figure 1 -** Participant Characteristics (based on
[^2^](https://www.zotero.org/google-docs/?R2k9SR))

Following additional quality control measures, raw metagenomic taxa were
quantified using MetaPhlAn4. To achieve this, MetaPhlAn4 compares the
acquired reads with a curated database consisting of 1.01 million
prokaryotic reference and metagenome-assembled genomes. The reads that
match these marker genes are then identified and quantified within the
sample, and each fragment is assigned a taxonomic label according to the
match [^12^](https://www.zotero.org/google-docs/?qXY1tC).

Further downstream analysis and visualization was conducted in R
(version 4.3.1), using the vegan, ecodist, fossil, and ggplot2 packages.

**2.2 Statistical Analyses**

**2.2.1 Diversity Analysis**

Alpha and beta diversity analyses were conducted to evaluate diversity
within each sample and group. In alpha-diversity analysis, relative
abundance values of species were evaluated using the Shannon, Simpson,
and Chao1 metrics, providing insights into species richness, evenness,
and estimated total richness. For beta-diversity analysis, Bray-Curtis
dissimilarity scores were used to quantify the dissimilarities between
control and disease patients, and these results were visualized using
Principal Coordinates Analysis (PCoA).

**2.2.2 Differential Abundance Analysis**

Differential abundance analysis across different groups was conducted
firstly by evaluating relative abundance at both species and genus
levels to identify the key contributors to liver disease. The Wilcoxon
rank-sum test was applied to these taxa, to subset those with p \< 0.05.
To further refine the subset, mean relative abundance values for each
taxon were calculated to display those that had the highest values and
categorize the rest as 'Others'. It is important to acknowledge the
trade-off of this approach - focusing on highly abundant taxa provides a
comprehensive view of significant contributors to liver disease but
sacrifices information on less abundant key players.

Then, the differences in mean relative abundances for species with p \<
0.05 were computed to provide insight into species-level changes
observed in liver disease patients versus control patients. The top 20
species with the most significant absolute changes in mean relative
abundance were then displayed, along with their means, and key species
were identified.

Finally, relative abundance analyses for the identified key species were
performed to confirm their significance. Differences across different
health states were observed, and the Wilcoxon rank-sum test was applied
to assess the statistical significance of these variations.

**2.2.3 Analysis of Confounding Variables**

To assess the impact of associated factors in liver disease, analyses of
age, gender, body mass index (BMI), albumin (Alb), creatinine (Crea),
total bilirubin (TB), and hepatitis B virus (HBV) were performed using
alpha diversity & Model for End-Stage Liver Disease (MELD) scores.

Shannon and Simpson metrics were employed to calculate alpha-diversity
for each factor, and Spearman's correlation coefficients were used to
determine the relations between the listed factors, diversity, and MELD
scores.

**3 Results**

**3.1 Diversity of the Gut Microbiome in Liver Disease**

Gut microbiome diversity for each sample was assessed through
alpha-diversity, employing the Simpson, Shannon, and Chao1 metrics. The
analysis of these findings unveils a distinct decrease in microbial
diversity within the gut of the diseased group across various metrics
(Figure 2).

Figure 2 provides valuable insights into the precision of each metric in
calculating alpha diversity for this dataset. Upon reviewing the
p-values calculated by the Wilcoxon rank-sum test for each metric, it
can be concluded that both the Simpson and Shannon metrics are suitable
for this dataset. The Simpson metric accounts for the dominance of a few
abundant species, while the Shannon metric considers both the number and
evenness of species. The Chao1 metric, however, exhibits a significantly
high p-value, especially in comparison to the other metrics. The p-value
associated with the Chao1 index substantially exceeds the conventional
threshold of p = 0.05 used to assess statistical significance,
suggesting that Chao1 may not be a sufficient metric for accurately
characterizing the alpha-diversity in this dataset.

![](media/image2.png){width="6.479034339457568in"
height="3.741058617672791in"}

**Figure 2 -** Boxplot representation of the alpha diversity of the gut
microbiome across control and disease groups, using Simpson, Shannon,
and Chao1 indexes (and Inverse Simpson for visual convenience). P-values
were calculated using the Wilcoxon rank-sum test.

The assessment of gut microbial diversity across groups was implemented
using the Bray-Curtis dissimilarity score and PCoA. Figure 3 showcases
the results obtained from these calculations, and provides insights into
how microbial community structures differ with the absence and presence
of liver disease.![](media/image3.png){width="4.618765310586177in"
height="3.843432852143482in"}

**Figure 3** - PcoA plot of Bray-Curtis dissimilarity of control and
liver disease groups with 95% confidence ellipses. Each dot represents a
measure of the microbiome composition of a given sample.

The tight clustering that can be observed in the healthy samples
suggests a similarity in microbial compositions within this group. In
contrast, the dispersion of points representing diseased samples
indicates higher variability in their microbial composition. The
increased microbial variability observed in the diseased samples may
suggest the presence of gut dysbiosis and could offer insights into the
severity of the disease.

**3.2 Identification of Key Species in Liver Disease**

Differential abundance analysis was conducted on the genus and species
data to pinpoint key species associated with the occurrence and
progression of liver disease. In Figure 4, the taxa with the highest
mean relative abundances are displayed at a genus and species level.
While genus-level analysis facilitates the identification of the major
groups present in the community and offers a broader perspective,
species-level analysis offers a more detailed view and aids the
identification of some potential key players in liver disease.

A\) B)

![](media/image4.png){width="3.057292213473316in"
height="3.057292213473316in"}![](media/image5.png){width="3.4739588801399823in"
height="2.989608486439195in"}

**Figure 4** - A) Stacked barplots representing the mean relative
abundance of genera present in liver disease and control patients B)
Stacked barplots representing the mean relative abundance of species
present in liver disease and control patients

Examining the genus-level data reveals a notable rise in the
*Veillonella* and *Streptococcus* genera, accompanied by a decline in
the *Alistipes* genus among patients with the disease. Similarly, at the
species level, we observe an upsurge in *Veillonella parvula*,
*Veillonella atypica*, and *Streptococcus salivarius*, coupled with a
reduction in *Alistipes putredinis* within the diseased group.
Additionally, though not apparent in its genus counterpart, a
significant decrease in *Bacteroides uniformis* is evident in this plot.

The species demonstrating the most significant shifts in mean relative
abundance during disease progression were graphically represented, along
with the mean differences, aiming for a more precise identification of
the pivotal species in liver disease (Figure 5). By cross-referencing
this information with the relative abundance figure presented above, we
can narrow down the list of key species to include *Alistipes
putredinis*, *Bacteroides uniformis*, *Veillonella parvula*,
*Streptococcus salivarius*, and *Veillonella atypica*.

When plotted, the relative abundance of these species reveals their
significant association with liver disease, as evidenced by the
corresponding p-values derived from the Wilcoxon rank-sum test (Figure
6).

![](media/image6.png){width="5.984375546806649in"
height="3.550729440069991in"}

**Figure 5** - Barplots illustrating species with the highest
differences in mean relative abundance, accompanied by the corresponding
mean difference values calculated as disease minus control.

![](media/image7.png){width="4.515625546806649in"
height="3.9979101049868766in"}

**Figure 6** - Boxplots showcasing the differential abundance of key
species in liver disease in control (red) and diseased patients (blue).

**3.3 Analysis of Factors Associated with Liver Disease**

Additional factors potentially influencing the gut microbiome were
evaluated in relation to liver disease. To achieve this, alpha diversity
correlation calculations were conducted for age, gender, BMI, Alb, Crea,
TB, and HBV data.

![](media/image8.png){width="5.395833333333333in"
height="3.350984251968504in"}

A\)

![](media/image9.png){width="5.3943110236220475in"
height="3.2604166666666665in"}

B\)

**Figure 7** - Scatter plots of alpha-diversity and factor correlation
using the Shannon and Inverse Simpson metrics. Each point represents a
healthy (red) or diseased (blue) sample. Regression lines and p-values
were calculated using Spearman's correlation coefficients.

The scatter plots depicting alpha-diversity vs. age, utilizing the
Shannon and Simpson metrics, reveal a negative slope for control
patients (Figure 7A). This indicates that, in general, an increase in
age is associated with decreased alpha-diversity. However, the barely
positive slope observed for diseased patients suggests a slight increase
in diversity with age. This trend remains inconclusive, and it may be
attributed to the complex interplay between the progression of liver
disease, the microbiota, and aging. Changes in liver function and the
immune system over time could potentially affect individuals with liver
disease differently as they age.

The negative slopes observed in both control and disease patients
underscore a distinct relationship between BMI and alpha-diversity
(Figure 7A). In both groups, higher BMI values are associated with
decreased alpha-diversity, a pattern supported by the low p-values
calculated with Spearman's coefficients.

While there doesn\'t appear to be significance with creatinine and total
bilirubin, the positive slope and low p-value observed in the alpha
diversity versus albumin plot indicate the increase and significance of
albumin levels in the context of liver disease (Figure 7B).

![](media/image10.png){width="2.6770833333333335in"
height="2.21662510936133in"}

A\) ![](media/image11.png){width="2.82874343832021in"
height="2.15625in"}

B\) C)![](media/image12.png){width="2.7707808398950133in"
height="2.15625in"}![](media/image13.png){width="2.6770833333333335in"
height="2.03125in"}

**Figure 8** - Boxplots illustrating the alpha-diversity of the gut
microbiome across A) female (red) and male (blue) groups B) HBV (red)
and no HBV (blue) groups C) alcohol (red) and no alcohol (blue) groups
according to Shannon and Inverse Simpson indices.

Gender, HBV, and alcohol consumption were also investigated in relation
to liver disease, but none of these factors stood out as significant, as
indicated by their high p-values (Figure 8).

In addition to alpha diversity correlation, MELD score correlation was
employed as another metric to assess the significance of the same
factors in liver disease. Upon examination of these scatter plots and
considering the insights gained from the alpha diversity plots, it can
be asserted that BMI and albumin exhibit a clear significance in the
occurrence and progression of liver disease.

![](media/image14.png){width="5.712889326334208in"
height="3.7551662292213472in"}

**Figure 9** - Scatter plots of age, BMI, creatinine, albumin, and TB in
correlation with MELD scores from diseased patients.

**4 Discussion**

This study utilized whole-metagenome shotgun sequencing to profile the
gut microbiome composition in both healthy individuals and those
afflicted with liver disease, aiming to identify key species and factors
associated with the condition. To achieve this, various statistical
analyses were performed including alpha-diversity analysis,
beta-diversity analysis, differential abundance analysis, and analysis
of confounding variables.

**4.1 Key Species in Liver Disease**

Findings from these analyses suggest elevated levels of *Veillonella
atypica*, *Veillonella parvula*, and *Streptococcus salivarius*, coupled
with a decrease in *Alistipes putredinis* and *Bacteroides uniformis*,
in the context of liver cirrhosis. These species are thereby emphasized
as key components in the dysbiotic gut microbiome of individuals with
liver disease.

General observations about the genera can also be made, suggesting
potential associations between liver cirrhosis and increased levels of
*Veillonella* and *Streptococcus*, along with a reduced abundance of the
*Alistipes* genus.

*Alistipes* is a recently identified genus known for its ability to
break down complex carbohydrates in the gut. Existing literature
suggests a connection between the *Alistipes* genus and protection
against certain diseases, including liver
fibrosis[^13^](https://www.zotero.org/google-docs/?QtA0AD). Notably, a
study by Iebba et. al (2018) indicated that a decrease in *Alistipes
spp*. correlates with the progression of liver cirrhosis into the
decompensated state [^14^](https://www.zotero.org/google-docs/?dwF4uW).
In liver cirrhosis, the levels of *A. shahii* and *A. putredinis* were
found to be decreased compared to healthy controls. Additionally, a
reduction in *Alistipes* abundance has been observed in patients with
liver fibrosis in other fibrotic diseases like NASH and
NAFLD[^15^](https://www.zotero.org/google-docs/?U5OYEI). Despite its
positive role in healthy phenotypes, *Alistipes* has been shown to have
a contrasting pathogenic role in diseases such as anxiety, myalgic
encephalomyelitis/chronic fatigue syndrome, and depression
[^13^](https://www.zotero.org/google-docs/?aZWIZe).

*Bacteroides uniformis* is typically a dominant species in the gut
microbiome, playing a crucial role in maintaining gut homeostasis, and
is abundantly present in healthy
individuals[^16,17^](https://www.zotero.org/google-docs/?1uQcz3). In
instances of liver disease, there is evidence suggesting a decrease in
the abundance of *B. uniformis*, contributing to the dysbiotic gut
microbiome observed in affected
patients[^16^](https://www.zotero.org/google-docs/?p45cPi). The low
abundance in liver disease patients can thus be attributed to the
dysbiotic composition in the gut, as evident by the alpha and beta
diversity analyses.

The increased abundance of *Streptococcus* and *Veillonella* points to a
particular trend in the translocation of oral bacteria to the gut, given
their prevalence in the oral microbiome of healthy patients
[^18^](https://www.zotero.org/google-docs/?SFZlsm). Recent studies
provide additional support for the translocation theory by indicating an
absence of evidence for the colonization of the gut by oral bacteria in
healthy individuals
[^19,20^](https://www.zotero.org/google-docs/?dUZeBG). The oral and gut
microbiomes are anatomically linked through saliva and food, acting as
mediums for this translocation. Kageyama et al. report aging as a
significant factor in the translocation of oral bacteria to the gut,
which prompts further exploration in the context of liver cirrhosis
given the inconclusive results obtained from the correlation of age and
liver disease [^20^](https://www.zotero.org/google-docs/?GUqZTe).

**4.2 Factors Associated with Liver Disease**

Typically, patients with liver cirrhosis exhibit reduced levels of
albumin, and human serum albumin (HSA) is a frequent method of treatment
for patients with liver cirrhosis
[^21^](https://www.zotero.org/google-docs/?BCRtxo). However, results
obtained from alpha diversity analysis (using Shannon and Inverse
Simpson metrics) and MELD score correlation analysis contradict this
phenomenon. Diseased patients in the data collected by Qin et al.
exhibit elevated albumin levels in contrast to their healthy
counterparts, indicating an additional complex interplay at play. Other
health conditions or medications may be influencing albumin levels
independently of liver function. Contribution to albumin production as
an inflammatory response to chronic inflammations is another possibility
for higher albumin levels in diseased patients
[^22^](https://www.zotero.org/google-docs/?pD6cIz).

Alpha diversity analysis and the correlation analysis with the Model for
End-Stage Liver Disease (MELD) score, conducted for Body Mass Index
(BMI), produce results consistent with existing literature. Elevated
BMIs are indicative of an elevated risk of liver disease, particularly
non-alcoholic fatty liver disease (NAFLD). A 2021 cross-sectional study
focusing on obese adolescents (BMI \> 30 kg/m2) concludes that
individuals with NAFLD and higher BMIs face an elevated risk of
developing liver fibrosis compared to those with lower BMIs
[^23^](https://www.zotero.org/google-docs/?L3EOmD). A meta-analysis
examining obesity as a prognostic factor for liver disease echoes
similar findings, highlighting the association between obesity and an
increased risk of severe liver
disease[^24^](https://www.zotero.org/google-docs/?pDt565).

**4.3 Future Works**

When navigating the trajectory for future research, it becomes crucial
to pinpoint gaps in knowledge, emerging trends, and evolving challenges
within the current research landscape. This section delves into
potential trajectories and areas of exploration that require further
investigation.

Future research should prioritize the sampling of whole-metagenomic
shotgun data with enhanced population diversity, as all the samples
collected by Qin et al. are exclusively from the Chinese population.
Increased population diversity offers several benefits, including a
broader understanding of the microbial landscape and the potential to
uncover variations across diverse genetic backgrounds and environmental
exposures.

Furthermore, further investigations aimed at comprehending the
biological mechanisms underlying the observed correlations could
elucidate the trends identified in the results. Conducting additional
investigations to unravel the biological mechanisms behind observed
correlations would provide insights into the involvement of immune and
metabolic factors and their influence on the correlated species.

And lastly, given the absence of an enhanced cure for liver cirrhosis,
it is crucial to contemplate alternative therapeutic approaches that can
inform and guide further research. Exploring strategies for targeting
the gut microbiome could offer valuable insights into preventing and
slowing disease progression.

**6 Conclusion**

In conclusion, this study sheds light on the relationship between the
gut microbiome and liver cirrhosis, identifying key species and factors
associated with the disease. Notably, the *Alistipes*, *Veillonella*,
*Streptococcus* genera, and *Bacteroides uniformis* play crucial roles
in the dysbiotic gut microbiome observed in liver disease patients.
Additionally, BMI and albumin emerged as significant factors,
challenging conventional expectations and demanding further exploration
for the latter. The findings emphasize the need for enhanced population
diversity in future research and underscore the importance of
investigating biological mechanisms and alternative therapeutic
approaches for a more holistic understanding and effective management of
liver cirrhosis.

**References**

[1. Lee, N. Y. & Suk, K. T. The Role of the Gut Microbiome in Liver
Cirrhosis Treatment. *Int. J. Mol. Sci.* **22**,
(2020).](https://www.zotero.org/google-docs/?KoCezR)

[2. Qin, N. *et al.* Alterations of the human gut microbiome in liver
cirrhosis. *Nature* **513**, 59--64
(2014).](https://www.zotero.org/google-docs/?KoCezR)

[3. Albillos, A., de Gottardi, A. & Rescigno, M. The gut-liver axis in
liver disease: Pathophysiological basis for therapy. *J. Hepatol.*
**72**, 558--577 (2019).](https://www.zotero.org/google-docs/?KoCezR)

[4. Smith, A., Baumgartner, K. & Bositis, C. Cirrhosis: Diagnosis and
Management. *Am. Fam. Physician* **100**, 759--770
(2019).](https://www.zotero.org/google-docs/?KoCezR)

[5. Jepsen, P. & Younossi, Z. M. The global burden of cirrhosis: A
review of disability-adjusted life-years lost and unmet needs. *J.
Hepatol.* **75**, S3--S13
(2021).](https://www.zotero.org/google-docs/?KoCezR)

[6. Guan, H., Zhang, X., Kuang, M. & Yu, J. The gut--liver axis in
immune remodeling of hepatic cirrhosis. *Front. Immunol.* **13**,
(2022).](https://www.zotero.org/google-docs/?KoCezR)

[7. Yan, A. W. *et al.* Enteric Dysbiosis Associated with a Mouse Model
of Alcoholic Liver Disease. *Hepathology* **53**, 96--105
(2011).](https://www.zotero.org/google-docs/?KoCezR)

[8. Bajaj, J., Heuman, D. & Hylemon, P. Altered profile of human gut
microbiome is associated with cirrhosis and its complications. *J
Hepatol* **940**, 940--947
(2013).](https://www.zotero.org/google-docs/?KoCezR)

[9. Bajaj, J. S. Altered Microbiota in Cirrhosis and Its Relationship to
the Development of Infection. *Clin. Liver Dis.* **14**, 107--111
(2019).](https://www.zotero.org/google-docs/?KoCezR)

[10. Loomba, R. *et al.* Gut Microbiome-Based Metagenomic Signature for
Non-invasive Detection of Advanced Fibrosis in Human Nonalcoholic Fatty
Liver Disease. *Cell Metab.* **25**, 1054--1062
(2017).](https://www.zotero.org/google-docs/?KoCezR)

[11. Durazzi, F. *et al.* Comparison between 16S rRNA and shotgun
sequencing data for the taxonomic characterization of the gut
microbiota. **11**, (2021).](https://www.zotero.org/google-docs/?KoCezR)

[12. Blanco-Miguez, A. *et al.* Extending and improving metagenomic
taxonomic profiling with uncharacterized species using MetaPhlAn4. *Nat.
Biotechnol.* (2023)
doi:https://doi.org/10.1038/s41587-023-01688-w.](https://www.zotero.org/google-docs/?KoCezR)

[13. Parker, B. J., Wearsch, P. A., Veloo, A. C. M. &
Rodriguez-Palacios, A. The Genus Alistipes: Gut Bacteria With Emerging
Implications to Inflammation, Cancer, and Mental Health. *Front.
Immunol.* **11**, (2020).](https://www.zotero.org/google-docs/?KoCezR)

[14. Iebba, V. *et al.* Combining amplicon sequencing and metabolomics
in cirrhotic patients highlights distinctive microbiota features
involved in bacterial translocation, systemic inflammation and hepatic
encephalopathy. *Nature*
doi:10.1038/s41598-018-26509-y.](https://www.zotero.org/google-docs/?KoCezR)

[15. Rau, M. *et al.* Fecal SCFAs and SCFA-producing bacteria in gut
microbiome of human NAFLD as a putative link to systemic T-cell
activation and advanced disease. *United Eur. Gastroenterol J.* **6**,
(2019).](https://www.zotero.org/google-docs/?KoCezR)

[16. Korobeinikova, A. V., Zlobovskaya, O. A. & Sheptulina, A. F. Gut
Microbiota Patterns in Patients with Non-Alcoholic Fatty Liver Disease:
A Comprehensive Assessment Using Three Analysis Methods. *Int. J. Mol.
Sci.* **24**, (2023).](https://www.zotero.org/google-docs/?KoCezR)

[17. Yan, Y. *et al.* Bacteroides uniformis-induced perturbations in
colonic microbiota and bile acid levels inhibit TH17 differentiation and
ameliorate colitis developments. *Nature* **9**,
(2023).](https://www.zotero.org/google-docs/?KoCezR)

[18. Yuan, H. *et al.* Quantitative changes of Veillonella,
Streptococcus, and Neisseria in the oral cavity of patients with
recurrent aphthous stomatitis: A systematic review and meta-analysis.
*Arch. Oral Biol.* **129**,
(2021).](https://www.zotero.org/google-docs/?KoCezR)

[19. Rashidi, A., Ebadi, M., Weisdorf, D. J., Costalonga, M. & Staley,
C. No evidence for colonization of oral bacteria in the distal gut in
healthy adults. *Proc. Natl. Acad. Sci. U. S. A.* **118**,
(2021).](https://www.zotero.org/google-docs/?KoCezR)

[20. Kageyama, S. *et al.* High-Resolution Detection of Translocation of
Oral Bacteria to the Gut. *J. Dent. Res.* **102**, 752--758
(2023).](https://www.zotero.org/google-docs/?KoCezR)

[21. Arroyo, V., Garcia-Martinez, R. & Salvatella, X. Human serum
albumin, systemic inflammation, and cirrhosis. *J. Hepatol.* **61**,
396--407 (2014).](https://www.zotero.org/google-docs/?KoCezR)

[22. Bernardi, M., Angeli, P., Claria, J. & Moreau, R. Albumin in
decompensated cirrhosis: new concepts and perspectives. *Gut*
**69**,.](https://www.zotero.org/google-docs/?KoCezR)

[23. Moran-Lev, H., Cohen, S., Webb, M. & Yerushalmy-Feler, A. Higher
BMI predicts liver fibrosis among obese children and adolescents with
NAFLD - an interventional pilot study. *BMC Pediatr.* **21**,
(2021).](https://www.zotero.org/google-docs/?KoCezR)

[24. Jarvis, H., Craig, D., Barker, R., Spiers, G. & Stow, D. Metabolic
risk factors and incident advanced liver disease in non-alcoholic fatty
liver disease (NAFLD): A systematic review and meta-analysis of
population-based observational studies. *PLOS Med.*
**17**,.](https://www.zotero.org/google-docs/?KoCezR)
