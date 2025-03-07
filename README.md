# FLa4a
- Version: 1.9.0.9002
- Author: Colin P. Millar and Ernesto Jardim
- Maintainer: Colin P. Millar <colin.millar AT ices.dk>
- Repository: <https://github.com/flr/FLa4a/>
- Bug reports: <https://github.com/flr/FLa4a/issues>

## Overview
FLa4a implements the Assesment For All (a4a) initiative stock assessment model, a simple and robust statistical catch-at-age model.

For further information about the a4a initiative, please visit <https://fishreg.jrc.ec.europa.eu/web/a4a>.

## Installation
To install this package, start R and enter:

    install.packages(c("copula","triangle", "coda"))

followed by

	install.packages("FLa4a", repos="http://flr-project.org/R")

or download from the [FLa4a releases page](https://github.com/flr/FLa4a/releases/latest)

## Documentation
- [Help pages](http://www.flr-project.org/FLa4a/reference/index.html)
- [Vignettes](http://www.flr-project.org/FLa4a/articles/index.html)

## References
- [Jardim, et.al, 2014](http://icesjms.oxfordjournals.org/content/early/2014/04/03/icesjms.fsu050.abstract)
- [Millar, et.al, 2014](http://icesjms.oxfordjournals.org/content/early/2014/03/31/icesjms.fsu043.abstract) 
- [Scott, et.al, 2016](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0154922)

## Build Status
[![Travis Build Status](https://travis-ci.org/flr/FLa4a.svg?branch=master)](https://travis-ci.org/flr/FLa4a)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/flr/FLa4a?branch=master&svg=true)](https://ci.appveyor.com/project/flr/FLa4)

## Releases
- [All releases](https://github.com/flr/FLCore/releases/)

## License
Copyright (c) 2012-2022 European Union. European Commission Joint Research Centre D.02. Released under the [EUPL 1.1](https://joinup.ec.europa.eu/community/eupl/home).

## Contact
You are welcome to:

- Submit suggestions and bug-reports at: <https://github.com/flr/FLa4a/issues>
- Send a pull request on: <https://github.com/flr/FLa4a/>
- Compose a friendly e-mail to the maintainer, see `packageDescription('FLa4a')`
