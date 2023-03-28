---
layout: default
title: Querying the SimSpin API
parent: Examples
nav_order: 4
last_modified_date: "Tues, 17 January 2022 11:11:00 AWST"
---

# Not a fan of R? Why not use the API...

In 2022, SimSpin was awarded development time with [Astronomy Data and Computing Services (ADACS)](https://adacs.org.au/who-we-are/) to build a web interface for working with the code. This application is a performant React Single Page App communicating asynchronously with a RESTful API to access the SimSpin package. 

[Launch the app!](https://simspin.datacentral.org.au/app/){: .btn .btn-purple }
{: .lh-tight }

Using the API, it is possible to run SimSpin from Python. In the example below, we take you through the process of building and downloading a mock observation.
{: .fs-5 .fw-300 .pb-2 }
---

To begin, we need to use the `requests` package for Python3. If you do not already have this pacakge installed, you can get it using the line below:

```Python
pip install requests
```

