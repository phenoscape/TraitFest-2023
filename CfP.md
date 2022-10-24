# Phenoscape TraitFest 2023	-- Call for Participation

## Synopsis

The [Phenoscape](https://phenoscape.org) team is hosting a 3.5 day hands-on workshop, TraitFest 2023, in January 2023 at [RENCI](https://renci.org), in Chapel Hill, North Carolina. If you are interested in generating, discovering, linking to, recombining, or computing with machine-interpretable trait data, this is the workshop for you!

The event will bring together a diverse group of people to collaboratively design and work on their research interests to take advantage of and promote reuse of Phenoscape’s online evolutionary data resources, tools, infrastructure, and services. The event is designed as a hands-on unconference-style workshop. Participants will break into subgroups to collaboratively tackle self-selected projects.

To apply to participate in the event, please fill out the [TraitFest 2023 Application for Participation](https://docs.google.com/forms/d/e/1FAIpQLSdWE5kOShjeeiMCo7DLsmFsjHIPboMpxF3vrutd_1kwWv52Xw/viewform) by the end of November 7, 2022. Funds to cover travel expenses are available but limited, as is space. We expect to notify applicants about acceptance starting one week after the application due date.

## Background: Computable trait data

Trait data that is amenable to computational data science, including computation-driven discovery, remains relatively new to science. Efficiently repurposing, integrating, and mining the vast stores of trait data has long been hampered by the limited amount of data accessible online in standard formats and by the challenges involved with enabling machines to compute with data that is largely recorded in natural language. A variety of advances has begun to address these challenges, including powerful knowledge representation technologies, shared domain ontologies, very fast machine reasoners, and various tools accelerating the construction of large databases (knowledge graphs) of ontology-linked phenotype data. The results emerging from these advances provide new opportunities for computation-driven data science with trait data. For biodiversity research, a unique resource of data science-enabled trait data is the Phenoscape knowledge base (KB), a datastore of vertebrate morphological traits linked to terms drawn from formal ontologies and thus represented in a structured and fully computable form.

## About the event

The Phenoscape TraitFest 2023 event aims to promote innovative applications adopting the capabilities conveyed by the data in the Phenoscape KB and its corresponding semantics-enabled tools, algorithms, and infrastructure. To do so, the event endeavors to engage, in a hands-on way, potential users of and contributors to this data, as well as developers of methods, such as in comparative phylogenetics and other fields who can use such data and algorithms. In particular, the event aims to include users whose research questions and/or datasets stand to benefit from being afforded the same computable semantics-based capabilities; whose taxonomic communities have already developed the necessary baseline infrastructure of supporting data and ontologies, such as the insects (e.g. Coleoptera, Hymenoptera, Diptera), spiders, and plant research communities; and those interested in developing tools or workflows, including in particular those using machine learning.

### Motivation

The overarching motivation for holding this event is to increase awareness of Phenoscape tools and resources, and to provide opportunities for hands-on knowledge exchange for using the data and computational capabilities, as well as applying the corresponding technologies and mechanisms to other taxonomic groups, trait data, and trait ontologies.

Particular potential outcomes we aim to facilitate include (but are not limited to) the following:
* Incubate new research or tool development projects, with enough momentum to continue after the event
* Learn about using and applying Phenoscape tools and technologies (such as Phenex, or Rphenoscape), and inform their continued development or improvement
* Inform and explore how Phenoscape tools and technologies can best be generalized to other taxonomic groups and sources of traits (e.g., images)
* Foster coordination between other anatomy and trait ontologies, and facilitate integration into the Phenoscape framework
* Explore integration or interoperability with other trait curation and annotation platforms, such as TaxonWorks
* Fostering new collaborations enabled by ontology-linked trait data, for example integrating morphological data into species distribution models, and integrating trait biology into environmental and other fields such as remote sensing

### Scope

We are keeping the scope of possible projects broad so as not to limit participants' ideas a priori. That said, we generally expect projects to use data and/or the computational services of Phenoscape's online resources (see “About the KB” below) in some way.

Project examples include, but are not limited to the following:

* Connecting one's data, tool, or database to the Phenoscape KB (to address a research question, or to interconnect a community resource)
* Increasing visibility into the KB data, such as through visualization, mashups, or apps
* Devising resources that make it easier to reuse the KB's data or capabilities, for example in the form of vignettes

Projects that align well with participants’ own professional interests are desirable as we hope these projects will be continued in some way after the event.

### Date and Location

The event will be held at the Renaissance Computing Institute (RENCI) in Chapel Hill, North Carolina. The exact date will be determined by polling participants’ availability. Current options are January 9-12 and January 23-26, 2023.

Updates will be posted to the workshop repository on Github: https://github.com/phenoscape/TraitFest-2023

## Who Should Participate

We are looking to assemble a diverse group of people, including domain experts from a variety of relevant fields, such as evolutionary biology, ecology, biodiversity science, biomedical sciences, bioinformatics, data science, machine learning and computer science. Examples of personas we are looking for include, but aren’t limited to the following:

* Researcher in a relevant field
* Software / tool / resource / data product developer
* Training / documentation specialist
* Visual and data interaction designer

Everyone participating in the event must agree to adhere to its [Code of Conduct](https://github.com/phenoscape/TraitFest-2023/blob/main/CODE_OF_CONDUCT.md).

## About Phenoscape, SCATE, and the Phenoscape KB

Since its inception in 2007, the [Phenoscape](https://phenoscape.org) project has worked to use ontologies and other knowledge representation and discovery technologies to render the natural language descriptions of comparative traits computable, and thus recombinable and interoperable. The project’s outcomes include resources such as the Phenoscape Knowledgebase (KB); tools for curating and staging data and ontologies for ingest; and visual and programmatic query interfaces. It’s current subproject, [SCATE](https://scate.phenoscape.org/) (Semantics for Comparative Analysis of Trait Evolution) is creating tools to enable KB data and computational capabilities assist in evolutionary analyses of trait evolution, in particular by providing comparative trait analysis tools easy access to algorithms powered by machine reasoning with the semantics of trait descriptions.

The [Phenoscape KB](http://kb.phenoscape.org) is an online resource containing evolutionary-relevant phenotypic trait data from more than 250 comparative morphology studies, with a current focus on the vertebrate fin-to-limb transition and comparative fish morphology. The KB contains natural language phenotype descriptions annotated with terms from formal ontologies (i.e., logical descriptions of phenotype concepts linked in a graph) so that machines can understand and compute on semantic (ontological) descriptions at scale. Using ontologies for morphological, spatial, and other requisite domain knowledge, the KB links these data to phenotypes reported in genetic perturbation studies of model organisms (zebrafish, frog, mouse) and to human genetic disease phenotypes. This enables discovery as the logic underlying ontologies can surface relationships between traits that were previously unknown using traditional methods.

The KB offers an application programming interface (API) to enable computational access to its trait datastore and exploration of connections among traits and between taxonomic groups. For querying, the KB uses machine reasoning to show both asserted (described in the literature) and inferred (connection determined by reasoning) traits of a taxonomic group. The API also gives access to other machine reasoning-based algorithms, such as inferring characters and states that are implied, but not necessarily asserted from original studies and finding evolutionary phenotype transitions semantically similar to gene phenotypes.
	
## Organizing Team

* Meghan Balk (Battelle, BGNN, Imageomics)
* Wasila Dahdul (UC Irvine, SCATE co-PI)
* Jennifer Girón Duque (Texas Tech University, Insect anatomy ontologies)
* Hilmar Lapp (Duke University, SCATE co-PI)
* Christopher Lawrence (Princeton University, Imageomics)
* István Mikó (University of New Hampshire, Insect anatomy ontologies)
