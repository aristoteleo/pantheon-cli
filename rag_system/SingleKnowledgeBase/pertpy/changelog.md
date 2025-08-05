# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.1


### 🚀 Features

* Add support for gamma parameter in MMD distance ([#825](https://github.com/scverse/pertpy/pull/825)) @Zethson
* Run notebooks in CI ([#815](https://github.com/scverse/pertpy/pull/815)) @Zethson

## 🐛 Bug Fixes

* Fix milo output writing ([#821](https://github.com/scverse/pertpy/pull/821)) @Zethson

## v1.0.0

### 🚀 Features

* Simplify guide assignment tutorial ([#796](https://github.com/scverse/pertpy/pull/796)) @Zethson
* Add milo pydeseq2 solver ([#795](https://github.com/scverse/pertpy/pull/795)) @Zethson

## v0.11.5

### 🚀 Features

* Move to trusted publishing ([#790](https://github.com/scverse/pertpy/pull/790)) @Zethson
* Add conda-forge to installation instructions ([#789](https://github.com/scverse/pertpy/pull/789)) @Zethson
* Fix seaborn 0.14 compat & testing dependencies & speed up tests ([#787](https://github.com/scverse/pertpy/pull/787)) @Zethson
* Support > 2 colors in `pt.tl.Sccoda.plot_boxplots` ([#784](https://github.com/scverse/pertpy/pull/784)) @mschilli87

### 🐛 Bug Fixes

* Fix LIL support for AnnData 0.12 ([#792](https://github.com/scverse/pertpy/pull/792)) @Zethson
* Fix compatibility with rpy 3.6.0 in milo ([#791](https://github.com/scverse/pertpy/pull/791)) @seohyonkim

## v0.11.4

### 🧰 Maintenance

* Remove decoupler ([#783](https://github.com/scverse/pertpy/pull/783)) @Zethson

## v0.11.3

### 🚀 Features

* add about page ([#770](https://github.com/scverse/pertpy/pull/770)) @Zethson
* Simplify Metadata errors ([#765](https://github.com/scverse/pertpy/pull/765)) @Zethson
* Standardize scCODA plot palette interface ([#773](https://github.com/scverse/pertpy/pull/773)) @mschilli87

## v0.11.2

### 🚀 Features

* Simplify Metadata errors ([#765](https://github.com/scverse/pertpy/pull/765)) @Zethson

## v0.11.1

### 🐛 Bug Fixes

* Restructure documentation & fix scCODA not requiring tcoda extra ([#762](https://github.com/scverse/pertpy/pull/762)) @Zethson

## v0.11.0

### 🚀 Features

* Improve CI test setup & speed up mixscape ([#757](https://github.com/scverse/pertpy/pull/757)) @Zethson
* Speed up mixscape ([#756](https://github.com/scverse/pertpy/pull/756)) @Zethson
* Optimize mixscape ([#735](https://github.com/scverse/pertpy/pull/735)) @Zethson
* Cleaner mixscape ([#746](https://github.com/scverse/pertpy/pull/746)) @eroell
* Speedup Mixscape's lda function ([#754](https://github.com/scverse/pertpy/pull/754)) @eroell
* Skip milotests properly & cache metadata ([#752](https://github.com/scverse/pertpy/pull/752)) @Zethson
* Sparse threshold max assignment ([#744](https://github.com/scverse/pertpy/pull/744)) @Zethson
* Speedup import times ([#743](https://github.com/scverse/pertpy/pull/743)) @stefanpeidli
* Add sparse threshold guide assignment & fix docstrings ([#737](https://github.com/scverse/pertpy/pull/737)) @Zethson
* More robust download ([#721](https://github.com/scverse/pertpy/pull/721)) @Zethson
* Support Python 3.13 ([#720](https://github.com/scverse/pertpy/pull/720)) @Zethson
* PRISM drug response metadata annotation ([#716](https://github.com/scverse/pertpy/pull/716)) @Lilly-May
* Improve CellLine metadata module UX ([#717](https://github.com/scverse/pertpy/pull/717)) @JiriBruthans
* Switch from ete3 to ete4 ([#618](https://github.com/scverse/pertpy/pull/618)) @dengzq1234
* Move pynndescent import out of for-loop ([#761](https://github.com/scverse/pertpy/pull/761)) @eroell
* Cleaner CinemaOT ([#745](https://github.com/scverse/pertpy/pull/745)) @eroell

### 🐛 Bug Fixes

* Change how guide assignment model saves params from uns to var ([#758](https://github.com/scverse/pertpy/pull/758)) @stefanpeidli
* A ton of doc improvements & fix mudata warnings & Augur interface ([#759](https://github.com/scverse/pertpy/pull/759)) @Zethson
* Change cell_wise_metric from "sqeuclidean" to "euclidean" for edist ([#733](https://github.com/scverse/pertpy/pull/733)) @stefanpeidli
* Check only for numeric in statsmodels ([#725](https://github.com/scverse/pertpy/pull/725)) @grst
* Remove leidenalg from requirements ([#729](https://github.com/scverse/pertpy/pull/729)) @Zethson
* Remove pyenchant ([#728](https://github.com/scverse/pertpy/pull/728)) @Zethson
* Fix coda plot axes ([#739](https://github.com/scverse/pertpy/pull/739)) @mschilli87

### 🧰 Maintenance

* A ton of doc improvements & fix mudata warnings & Augur interface ([#759](https://github.com/scverse/pertpy/pull/759)) @Zethson
* Additional pre-commit checks ([#748](https://github.com/scverse/pertpy/pull/748)) @Zethson
* Mixscape dropdown ([#740](https://github.com/scverse/pertpy/pull/740)) @Zethson
* Fix docstrings ([#737](https://github.com/scverse/pertpy/pull/737)) @Zethson

## v0.10.0

### 🚀 Features

* Allow Custom Feature Spaces in Dialogue ([#712](https://github.com/scverse/pertpy/pull/712)) @grpinto
* Add Poisson-Gauss Mixture model for Guide Assignment ([#709](https://github.com/scverse/pertpy/pull/709)) @stefanpeidli
* Add perturbation signature calculation from replicate control cells ([#695](https://github.com/scverse/pertpy/pull/695)) @Lilly-May
* Optimize test speed ([#699](https://github.com/scverse/pertpy/pull/699)) @Zethson
* Optimize MeanVarDistributionDistance ([#697](https://github.com/scverse/pertpy/pull/697)) @Zethson

### 🐛 Bug Fixes

* Align Mixscape with Seurat's implementation ([#710](https://github.com/scverse/pertpy/pull/710)) @Lilly-May
* Fix empty figure returns when show=True in plotting functions ([#703](https://github.com/scverse/pertpy/pull/703)) @Lilly-May
* Fix probability data type ([#696](https://github.com/scverse/pertpy/pull/696)) @Zethson
* Mixscape reproducibility ([#683](https://github.com/scverse/pertpy/pull/683)) @Lilly-May

## v0.9.x

### 🚀 Features

* Python 3.12 support ([#644](https://github.com/scverse/pertpy/pull/644)) @Lilly-May
* Add lazy loading for deferred modules ([#647](https://github.com/scverse/pertpy/pull/647) and [#648](https://github.com/scverse/pertpy/pull/648)) @Zethson
* Support for scipy 1.14 ([#646](https://github.com/scverse/pertpy/pull/646)) @andr-kun
* Update scperturb datasets ([#649](https://github.com/scverse/pertpy/pull/649)) @tessadgreen
* Augur scsim warnings ([#670](https://github.com/scverse/pertpy/pull/670)) @Zethson
* Add uncertainty score in KNN label_transfer in PerturbationSpace ([#658](https://github.com/scverse/pertpy/pull/658)) @stefanpeidli
* Incorporate use case tutorials ([#665](https://github.com/scverse/pertpy/pull/665)) @Lilly-May
* Add plots for DE analysis and unify plotting API ([#654](https://github.com/scverse/pertpy/pull/654)) @Lilly-May
* Switch to formulaic-contrasts ([#682](https://github.com/scverse/pertpy/pull/682)) @grst @emdann
* Remove anndata pin ([#653](https://github.com/scverse/pertpy/pull/653)) @Lilly-May

### 🐛 Bug Fixes

* Fix jax random array ([#686](https://github.com/scverse/pertpy/pull/686)) @Zethson
* Fixed plotting for mixscape.plot_barplot and sccoda.plot_effects_barplot ([#667](https://github.com/scverse/pertpy/pull/667)) @Lilly-May
* Remove explicit anndata in dependencies ([#666](https://github.com/scverse/pertpy/pull/666)) @Zethson
* Fix legend position for pt.tl.Mixscape.plot_barplot ([#660](https://github.com/scverse/pertpy/pull/660)) @Lilly-May

### 🧰 Maintenance

* Fix docs rendering for classes using lazy import ([#651](https://github.com/scverse/pertpy/pull/651)) @Lilly-May

## v0.8.0

### 🚀 Features

* Add hagai dataloader ([#627](https://github.com/scverse/pertpy/pull/627)) @Zethson
* Add more informative error if not every cell type is represented in e… ([#620](https://github.com/scverse/pertpy/pull/620)) @Zethson
* Add typehints_defaults ([#619](https://github.com/scverse/pertpy/pull/619)) @Zethson
* Add metadata notebook ([#609](https://github.com/scverse/pertpy/pull/609)) @Zethson
* Add multicondition analysis (differential gene expression interface) ([#607](https://github.com/scverse/pertpy/pull/607)) @Zethson
* Add zhang 2021 dataloader ([#602](https://github.com/scverse/pertpy/pull/602)) @Zethson
* Speed up tests ([#601](https://github.com/scverse/pertpy/pull/601) [#599](https://github.com/scverse/pertpy/pull/599)) @Zethson
* Use __all__ ([#588](https://github.com/scverse/pertpy/pull/588)) @Zethson
* Use lamin logger ([#586](https://github.com/scverse/pertpy/pull/586)) @Zethson
* Add two distance metrics, three-way comparison and bootstrapping ([#608](https://github.com/scverse/pertpy/pull/608)) @wxicu
* Metadata GDSC annotation improvements ([#625](https://github.com/scverse/pertpy/pull/625)) @Lilly-May
* Add file lock during downloading ([#590](https://github.com/scverse/pertpy/pull/590)) @wxicu

### 🐛 Bug Fixes

* Fix CI ([#637](https://github.com/scverse/pertpy/pull/637)) @Zethson
* Fix layer handling in Mixscape.mixscape to resolve errors with adata.raw ([#636](https://github.com/scverse/pertpy/pull/636)) @Lilly-May
* Fix from_scanpy ([#596](https://github.com/scverse/pertpy/pull/596)) @Zethson
* Fix volcano plot y-axis being cut-off ([#622](https://github.com/scverse/pertpy/pull/622)) @pakiessling

## v0.7.0

Due to a lack of release notes for 0.6.0, this release may contain the some of the changes of 0.6.0 as well.

Note that the syntax for plotting for most tools changed. We removed the `pt.pl` modules in favor of moving all plots with `plot_` prefixes to the respective tools. We hope that this will make the documentation easier to navigate.

### 🚀 Features

* Silence pandas groupby() warning when running Augur predict ([#573](https://github.com/scverse/pertpy/pull/573)) @namsaraeva
* Added milo beeswarm example plot ([#552](https://github.com/scverse/pertpy/pull/552)) @namsaraeva
* Lazily metadata loading ([#544](https://github.com/scverse/pertpy/pull/544)) @wxicu
* Added tasccoda example plots ([#527](https://github.com/scverse/pertpy/pull/527)) @namsaraeva
* Add examples to plotting functions ([#511](https://github.com/scverse/pertpy/pull/511)) @namsaraeva
* Replace plotnine implementation of barplot with seaborn/matplotlib ([#441](https://github.com/scverse/pertpy/pull/441)) @kcArtemis
* Cinemaot dataset loader ([#424](https://github.com/scverse/pertpy/pull/424)) @Lilly-May
* Add classifiers ([#427](https://github.com/scverse/pertpy/pull/427)) @Zethson
* Remove remotezip ([#401](https://github.com/scverse/pertpy/pull/401)) @wxicu
* Add new distances ([#304](https://github.com/scverse/pertpy/pull/304)) @mojtababahrami
* wasserstein distance return type float ([#386](https://github.com/scverse/pertpy/pull/386)) @eroell
* Implementation of CINEMA-OT for pertpy - follow up ([#379](https://github.com/scverse/pertpy/pull/379)) @MingzeDong
* Implementation of CINEMA-OT for pertpy ([#377](https://github.com/scverse/pertpy/pull/377)) @MingzeDong
* Update cell line meta data class ([#539](https://github.com/scverse/pertpy/pull/539)) @wxicu
* Add dgidb and pharmgkb for drug annotation ([#575](https://github.com/scverse/pertpy/pull/575)) @wxicu
* Add sciplex-GxE dataloader ([#576](https://github.com/scverse/pertpy/pull/576)) @namsaraeva
* Harmonize plots ([#574](https://github.com/scverse/pertpy/pull/574)) @namsaraeva
* Compute method for MLPClassifierSpace ([#565](https://github.com/scverse/pertpy/pull/565)) @Lilly-May
* uv for CI ([#569](https://github.com/scverse/pertpy/pull/569)) @Zethson
* Logistic regression support for the Discriminator Classifier ([#560](https://github.com/scverse/pertpy/pull/560)) @Lilly-May
* Allow saving of CODA plots ([#559](https://github.com/scverse/pertpy/pull/559)) @namsaraeva
* Order plot functions last ([#555](https://github.com/scverse/pertpy/pull/555)) @Zethson
* Sparse guide RNA plot ([#554](https://github.com/scverse/pertpy/pull/554)) @Zethson
* Rename cmap to palette ([#553](https://github.com/scverse/pertpy/pull/553)) @Zethson
* Expand load for sccoda docstrings ([#543](https://github.com/scverse/pertpy/pull/543)) @Zethson
* Pseudobulk improvements ([#529](https://github.com/scverse/pertpy/pull/529)) @Lilly-May
* Add more detailed tool explanations ([#528](https://github.com/scverse/pertpy/pull/528)) @Zethson
* Add knn imputation ([#517](https://github.com/scverse/pertpy/pull/517)) @Zethson
* Add combosciplex ([#512](https://github.com/scverse/pertpy/pull/512)) @Zethson
* Add DE set comparisons ([#489](https://github.com/scverse/pertpy/pull/489)) @Zethson
* More informative images in tutorial gallery ([#488](https://github.com/scverse/pertpy/pull/488)) @Lilly-May
* Add perturbation space tutorial ([#487](https://github.com/scverse/pertpy/pull/487)) @Lilly-May
* Add enrichment ([#482](https://github.com/scverse/pertpy/pull/482)) @Zethson
* Save perturbation and clustering labels as categorical ([#481](https://github.com/scverse/pertpy/pull/481)) @Lilly-May
* Optimize from_scanpy ([#473](https://github.com/scverse/pertpy/pull/473)) @Zethson
* Add Drug metadata from chembl ([#480](https://github.com/scverse/pertpy/pull/480)) @Zethson
* Add DE filtering ([#477](https://github.com/scverse/pertpy/pull/477)) @Zethson
* PertSpace docs improvements ([#471](https://github.com/scverse/pertpy/pull/471)) @Lilly-May
* Make dialogue's load protected ([#464](https://github.com/scverse/pertpy/pull/464)) @Lilly-May
* New plotting API ([#456](https://github.com/scverse/pertpy/pull/456)) @Zethson
* Added smillie_2019 dataloader and renamed smillie dataset to tasccoda_example ([#450](https://github.com/scverse/pertpy/pull/450)) @Lilly-May
* Docs examples for CINEMA-OT ([#433](https://github.com/scverse/pertpy/pull/433)) @Lilly-May
* Removing dependencies ([#426](https://github.com/scverse/pertpy/pull/426)) @Zethson
* Added example for pt.pl.milo.nhood_graph ([#423](https://github.com/scverse/pertpy/pull/423)) @Lilly-May
* scope imports of ete3 in tasccoda ([#422](https://github.com/scverse/pertpy/pull/422)) @Zethson
* Add cinemaot nb ([#419](https://github.com/scverse/pertpy/pull/419)) @Zethson
* DIALOGUE attribution ([#371](https://github.com/scverse/pertpy/pull/371)) @Zethson
* DIALOGUE extensions and plots ([#368](https://github.com/scverse/pertpy/pull/368)) @tessadgreen
* Set run nuts random key to 0 by default. ([#363](https://github.com/scverse/pertpy/pull/363)) @Zethson
* Remove statsannotations ([#362](https://github.com/scverse/pertpy/pull/362)) @Zethson
* Refactoring ([#354](https://github.com/scverse/pertpy/pull/354)) @Zethson
* Add combosciplex ([#512](https://github.com/scverse/pertpy/pull/512)) @Zethson

### 🐛 Bug Fixes

* Fix stephenson 2021 and refactor ([#374](https://github.com/scverse/pertpy/pull/374)) @Zethson
* fix: double plot issue related to barplot ([#460](https://github.com/scverse/pertpy/pull/460)) @kcArtemis
* Small fixes to milo.plot_da_beeswarm ([#551](https://github.com/scverse/pertpy/pull/551)) @emdann
* Check and fix plotting functions ([#579](https://github.com/scverse/pertpy/pull/579)) @namsaraeva
* Fix jax config ([#550](https://github.com/scverse/pertpy/pull/550)) @Zethson
* fixed augur runtime ([#547](https://github.com/scverse/pertpy/pull/547)) @zzheng18
* Fix Distance docstring rendering ([#537](https://github.com/scverse/pertpy/pull/537)) @Zethson
* Address some Jax random key issues ([#531](https://github.com/scverse/pertpy/pull/531)) @Zethson
* Fixed Augur bug ([#533](https://github.com/scverse/pertpy/pull/533)) @namsaraeva
* Fix discriminator classifier nn dimensions ([#475](https://github.com/scverse/pertpy/pull/475)) @Lilly-May
* bug fixes, implementation tweaks, and additional distances ([#397](https://github.com/scverse/pertpy/pull/397)) @yugeji
* Fix subsetting of milo.plot_da_beeswarm ([#472](https://github.com/scverse/pertpy/pull/472)) @Lilly-May
* CentroidSpace AnnData Annotations ([#455](https://github.com/scverse/pertpy/pull/455)) @Lilly-May
* DiscriminatorClassifier compatibility with sparse matrices ([#453](https://github.com/scverse/pertpy/pull/453)) @Lilly-May
* Fix code coverage ([#432](https://github.com/scverse/pertpy/pull/432)) @flying-sheep
* Fixed set_fdr for tascCODA models ([#411](https://github.com/scverse/pertpy/pull/411)) @Lilly-May
* Small bug fixes and docs improvements ([#409](https://github.com/scverse/pertpy/pull/409)) @Lilly-May
* Move to new numpy.random.Generator  ([#396](https://github.com/scverse/pertpy/pull/396)) @Lilly-May
* Fix cinema OT test ([#392](https://github.com/scverse/pertpy/pull/392)) @Zethson
* fix naming of example data in doc examples ([#387](https://github.com/scverse/pertpy/pull/387)) @eroell
* Fix stephenson 2021 and refactor ([#374](https://github.com/scverse/pertpy/pull/374)) @Zethson
* scCODA: add low/high acceptance probability warnings ([#366](https://github.com/scverse/pertpy/pull/366)) @johannesostner
* actually use `neighbors_key` in `milo` ([#418](https://github.com/scverse/pertpy/pull/418)) @maarten-devries
* Fixed some bugs (due to inconsistencies during migration) and typos in CINEMA-OT ([#413](https://github.com/scverse/pertpy/pull/413)) @MingzeDong
* fix Implicit modification warnings in cellline tests ([#404](https://github.com/scverse/pertpy/pull/404)) @wxicu

## v0.5.0

### 🚀 Features

* Perturbation space - add/subtract ([#328](https://github.com/scverse/pertpy/pull/328)) @AlejandroTL
* Perturbation space - pseudobulk, discriminative classifier, clustering ([#316](https://github.com/scverse/pertpy/pull/316)) @AlejandroTL
* add n_jobs support for pairwise distance computation ([#305](https://github.com/scverse/pertpy/pull/305)) @XinmingTu

### 🐛 Bug Fixes

* fix initial MCMC values ([#302](https://github.com/scverse/pertpy/pull/302)) @johannesostner
* Fix docs CI ([#303](https://github.com/scverse/pertpy/pull/303)) @Zethson
* fixed the output_file_name of frangieh_2021_raw ([#321](https://github.com/scverse/pertpy/pull/321)) @XinmingTu

## v0.4.0

### 🚀 Features

* Guide assignment ([#194](https://github.com/scverse/pertpy/pull/194)) @moinfar
* Distances ([#217](https://github.com/scverse/pertpy/pull/217)) @Zethson
* Add initial DIALOGUE ([#199](https://github.com/scverse/pertpy/pull/199)) @Zethson
* Add Jax ScGen ([#183](https://github.com/scverse/pertpy/pull/183)) @Zethson
* Guide assignment ([#194](https://github.com/scverse/pertpy/pull/194)) @moinfar
* Add cell line annotation ([#187](https://github.com/scverse/pertpy/pull/187)) @wxicu

### 🐛 Bug Fixes

* docstring enhancement for sample identifiers ([#237](https://github.com/scverse/pertpy/pull/237)) @johannesostner
* improve mixscape documentation & fix a bug ([#222](https://github.com/scverse/pertpy/pull/222)) @xinyuejohn
* Fix docs ete3 ([#206](https://github.com/scverse/pertpy/pull/206)) @Zethson
* Update _dialogue.py ([#203](https://github.com/scverse/pertpy/pull/203)) @tessadgreen
* scCODA bugfix: duplicates in covariate_obs ([#198](https://github.com/scverse/pertpy/pull/198)) @Zethson
* Fix [#211](https://github.com/scverse/pertpy/issues/211) ([#212](https://github.com/scverse/pertpy/pull/212)) @Moomboh
* fix [#207](https://github.com/scverse/pertpy/issues/207) ete3 warning masks other ImportErrors ([#208](https://github.com/scverse/pertpy/pull/208)) @Moomboh

## v0.3.0

### 🚀 Features

* added list option for sample_identifier in load ([#190](https://github.com/scverse/pertpy/pull/190)) @Zethson
* integrate scCODA and tascCODA ([#188](https://github.com/scverse/pertpy/pull/188)) @xinyuejohn
* add milopy ([#165](https://github.com/scverse/pertpy/pull/165)) @xinyuejohn

### 🐛 Bug Fixes

* changes default settings of coda and milo ([#196](https://github.com/scverse/pertpy/pull/196)) @xinyuejohn
* fix mixscape default setting of layer ([#195](https://github.com/scverse/pertpy/pull/195)) @xinyuejohn

### 🧰 Maintenance

* integrate scCODA and tascCODA ([#188](https://github.com/scverse/pertpy/pull/188)) @xinyuejohn
* add milopy ([#165](https://github.com/scverse/pertpy/pull/165)) @xinyuejohn

## v0.2.0

### 🚀 Features

* added list option for sample_identifier in load ([#190](https://github.com/scverse/pertpy/pull/190)) @Zethson
* integrate scCODA and tascCODA ([#188](https://github.com/scverse/pertpy/pull/188)) @xinyuejohn
* add milopy ([#165](https://github.com/scverse/pertpy/pull/165)) @xinyuejohn
* Add scperturb datasets ([#182](https://github.com/scverse/pertpy/pull/182)) @wxicu
* Added subsampled Stephenson2021 COVID19 dataset ([#179](https://github.com/scverse/pertpy/pull/179)) @emdann

### 🐛 Bug Fixes

* changes default settings of coda and milo ([#196](https://github.com/scverse/pertpy/pull/196)) @xinyuejohn
* fix mixscape default setting of layer ([#195](https://github.com/scverse/pertpy/pull/195)) @xinyuejohn
* Solve pre-commit prettier and isort version incompatibility issues ([#193](https://github.com/scverse/pertpy/pull/193)) @moinfar

### 🧰 Maintenance

* integrate scCODA and tascCODA ([#188](https://github.com/scverse/pertpy/pull/188)) @xinyuejohn
* add milopy ([#165](https://github.com/scverse/pertpy/pull/165)) @xinyuejohn

## v0.1.0

### 🚀 Features

* Augurpy refactoring ([#148](https://github.com/scverse/pertpy/pull/148)) @Zethson
* add kang dataset ([#146](https://github.com/scverse/pertpy/pull/146)) @antschum
* add inits for plots, fix error messages and output, tutorial version ([#144](https://github.com/scverse/pertpy/pull/144)) @antschum
* Feature/furo md ([#141](https://github.com/scverse/pertpy/pull/141)) @Zethson
* add augurpy ([#135](https://github.com/scverse/pertpy/pull/135)) @antschum
* add better data draft ([#71](https://github.com/scverse/pertpy/pull/71)) @Zethson
* Documentation structure ([#23](https://github.com/scverse/pertpy/pull/23)) @Zethson
* add pypi-latest & architecture ([#22](https://github.com/scverse/pertpy/pull/22)) @Zethson

### 🐛 Bug Fixes

* add inits for plots, fix error messages and output, tutorial version ([#144](https://github.com/scverse/pertpy/pull/144)) @antschum

### 🧰 Maintenance

* Augurpy refactoring ([#148](https://github.com/scverse/pertpy/pull/148)) @Zethson
* Feature/furo md ([#141](https://github.com/scverse/pertpy/pull/141)) @Zethson
* add pypi-latest & architecture ([#22](https://github.com/scverse/pertpy/pull/22)) @Zethson
