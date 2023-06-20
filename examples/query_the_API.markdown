---
layout: default
title: Querying the SimSpin API
parent: Examples
nav_order: 4
last_modified_date: "Tue, 20 June 2023 11:11:00 AWST"
---

# Not a fan of R? Why not use the API...

In 2022, SimSpin was awarded development time with [Astronomy Data and Computing Services (ADACS)](https://adacs.org.au/who-we-are/) to build a web interface for working with the code. This application is a performant React Single Page App communicating asynchronously with a RESTful API to access the SimSpin package. 

[Launch the app!](https://simspin.datacentral.org.au/app/){: .btn .btn-purple }
{: .lh-tight }

Using the API, it is possible to run SimSpin from Python. In the example below, we take you through the process of building and downloading a mock observation.
{: .fs-5 .fw-300 .pb-2 }
---

## Interacting with the API using `requests`

To begin, we need to use the `requests` package for Python3. If you do not already have this package installed, you can get it using the line below:

```python
pip install requests
```

With this package installed, we can go about working with the SimSpin package via queries to the API.

---

## Building a SimSpin file through the API

Let's begin by building a SimSpin file. Using an example file from the package, we can post a query to the `make_simspin_file` function using the following code snippet. Here we do not describe the input parameters for each function, but you can find out more about each one at the [`make_simspin_file`](/SimSpin/docs/make_simspin_file.markdown) documentation page. 

```python
import requests

# Declare the input parameters for running `make_simpsin_file`
parameters = {
    "email": "simspin@fastmail.com",
    "disk_age": 5,
    "bulge_age": 5,
    "disk_Z": 0.004,
    "bulge_Z": 0.004,
    "template": "BC03lr",
    "centre": 0,
    "half_mass": 0,
    "sph_spawn_n": 10
}

file = {"snapshot_file": open("SimSpin_example_EAGLE.hdf5", "rb")}

# Post the `make_simpsin_file` request to the API
post_SimSpin_file_request = requests.post(
    "https://simspin.datacentral.org.au/api/make_simspin_file/",
    data  = parameters,
    files = file
)
```
We can check on the resulting output of the function using:

```python
post_SimSpin_file_request.json()

{'id': 'a83ca9f9-823a-4739-9721-23fceacecdd8',
 'url': 'https://simspin.datacentral.org.au/api/make_simspin_file/a83ca9f9-823a-4739-9721-23fceacecdd8/',
 'created': '2023-05-23T06:14:38.856828Z',
 'email': 'simspin@fastmail.com',
 'snapshot_file': 'https://simspin.datacentral.org.au/storage/public/make_simspin_file/a83ca9f9-823a-4739-9721-23fceacecdd8/uploads/SimSpin_example_Magneticum.hdf5',
 'disk_age': 5.0,
 'bulge_age': 5.0,
 'disk_Z': 0.004,
 'bulge_Z': 0.004,
 'template': 'BC03lr',
 'centre': 0.0,
 'half_mass': 0.0,
 'sph_spawn_n': 10.0}

 ```
 
 We will also receive an email to the listed address once the file construction is complete, as shown below. 
<img align="centre" src="assets/images/simspin_file_api_example_email.png" width="750" height="237" />

If you are a registered Data Central user, you can login to your account to browse a history of constructed files. In this scenario, we have not specififed our user credentials so must interact with the file using the `id` link listed. Navigating to the URL listed against `url` allows us to browse the details of the constructed file and download it directly if you wish to explore the file locally. 

Otherwise, we can proceed to interacting with the produced SimSpin file on the server, and use this to construct a data cube.

---

## Building a data cube using the API

