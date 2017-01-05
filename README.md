# FLa4a
- Version: 1.0.0
- Date: 2016-12-21
- Author: Colin P. Millar and Ernesto Jardim
- Maintainer: Ernesto Jardim <ernesto.jardim AT jrc.ec.europa.eu>
- Repository: <https://github.com/flr/FLa4a/>
- Bug reports: <https://github.com/flr/FLa4a/issues>

## Overview
FLa4a implements the Assesment For All (A4A) initiative stock assessment model, a simple and robust statistical catch-at-age model. For further information about the A4A initiative, please visit <https://fishreg.jrc.ec.europa.eu/web/a4a>.

To install this package, start R and enter:

	install.packages(c("copula", "mgcv", "triangle", "coda"))

followed by

	install.packages("FLa4a", repos="http://flr-project.org/R")

or download from the [FLa4a releases page](https://github.com/flr/FLa4a/releases/latest)

## Documentation
- [Help pages](http://flr-project.org/FLa4a)

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
Copyright (c) 2012-2016 European Union. European Commission Joint Research Centre D.02. Released under the [EUPL 1.1](https://joinup.ec.europa.eu/community/eupl/home).

## Contact
You are welcome to:

- Submit suggestions and bug-reports at: <https://github.com/flr/FLa4a/issues>
- Send a pull request on: <https://github.com/flr/FLa4a/>
- Compose a friendly e-mail to the maintainer, see `packageDescription('FLa4a')`
