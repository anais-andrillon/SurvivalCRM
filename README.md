The codes in this directory can be used to reproduce the results in the paper

# "Dose-finding design and benchmark for a right censored endpoint"
by
# Anaïs Andrillon, Sylvie Chevret, Shing M Lee, and Lucie Biard (2020)

In our paper, we developped the survival-Continual Reassessment Method (#Surv-CRM#) design building on the CRM dose-finding design, but using survival models for right-censored DLT endpoints allowing the outcomes to be delayed. To handle possible informative censoring, we also proposed the informative survival-CRM (#iSurv-CRM#), that extends the Surv-CRM by considering a competing-risk model. In addition, we developed a nonparametric benchmark approach for evaluation of dose-finding designs with right censored time-to-event endpoints. 

This work was motivated by the need for specific methods for dose-finding with novel anti-cancer agents, such as targeted therapies or immunotherapies, which often require prolonged observation windows and result in some patients who do not experience a DLT but have not yet reached the end of the window. Moreover, when the observation window is long relatively to the underlying disease process, the observation of toxicity may be precluded by trial discontinuation related to lack of efficacy of the drug (due to death, progression, withdrawal, etc.).
 
## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
