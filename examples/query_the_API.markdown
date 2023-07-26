---
layout: default
title: Querying the SimSpin API
parent: Examples
nav_order: 4
last_modified_date: "Wed, 26 July 2023 11:11:00 AWST"
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

Let's begin by building a SimSpin file. Using an example file from the package, we can post a query to the `make_simspin_file` function using the following code snippet. We do not describe the input parameters for each function here, but you can find out more about each one at the [`make_simspin_file`](/SimSpin/docs/make_simspin_file) documentation page. 

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

{'id': '5faf1a8c-3225-4368-a2b8-e5c02f0359b9',
 'url': 'https://simspin.datacentral.org.au/api/make_simspin_file/5faf1a8c-3225-4368-a2b8-e5c02f0359b9/',
 'created': '2023-07-26T08:14:12.713279Z',
 'email': 'simspin@fastmail.com',
 'snapshot_file': 'https://simspin.datacentral.org.au/storage/public/make_simspin_file/5faf1a8c-3225-4368-a2b8-e5c02f0359b9/uploads/SimSpin_example_EAGLE.hdf5',
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

<img align="centre" src="assets/images/simspin_file_api_example_email.png" width="500" height="158" />

If you are a registered Data Central user, you can login to your account to browse a history of constructed files. In this scenario, we have not specififed our user credentials so must interact with the file using the `id` link listed. Navigating to the URL listed against `url` allows us to browse the details of the constructed file and download it directly if you wish to explore the file locally. 

Navigating to the [url link](https://simspin.datacentral.org.au/api/make_simspin_file/5faf1a8c-3225-4368-a2b8-e5c02f0359b9/) above, we see the following information about the created file and the `simspin_file` link to the public simspin file created.

<img align="centre" src="assets/images/simspin_file_api_example_web_interface.png" width="750" height="237" />

If you are following along with this example and would like to download a copy of the SimSpin file created, you can do so using the button below. Otherwise, we can proceed to interacting with the produced SimSpin file on the server, and use this to construct a data cube.

<span class="fs-4">
[Get the SimSpin file](https://simspin.datacentral.org.au/storage/public/make_simspin_file/5faf1a8c-3225-4368-a2b8-e5c02f0359b9/files/SimSpin_example_EAGLE_BC03lr.Rdata){: .btn .btn-purple }
</span>

---

## Building a data cube using the API

Proceeding with building the observed data cube using the SimSpin file constructed above, we organise the various inputs of the [`build_datacube`](/SimSpin/docs/build_datacube) function into a Python dictionary.

The code snippet below will download our simulated galaxy, project it to a distance of 0.3 kpc/'' at an inclination of 30&deg;and observe it with a SAMI-like instrument. 

We use the "." notation to specify the inputs for [`telescope`](/SimSpin/docs/telescope) and [`observing_strategy`](/SimSpin/docs/observing_strategy) within the data dictionary. 
{: .note}

```python
# Downloading and saving the SimSpin file
EAGLE_simspin_file = "https://simspin.datacentral.org.au/storage/public/make_simspin_file/5faf1a8c-3225-4368-a2b8-e5c02f0359b9/files/SimSpin_example_EAGLE_BC03lr.Rdata"
r = requests.get(EAGLE_simspin_file, allow_redirects=True)
open('SimSpin_example_EAGLE_BC03lr.Rdata', 'wb').write(r.content)

# Listing the inputs for the build_datacube function
input_files = {'simspin_file': open('SimSpin_example_EAGLE_BC03lr.Rdata', 'rb')}

data = {
    "method": "velocity",
    "object_name": "EAGLE_example_SAMI",
    "telescope_name": "SAMI",
    "observer_name": "SimSpin",
    "split_save": True,
    "mass_flag": False,
    "observing_strategy.dist_z": None,
    "observing_strategy.dist_Mpc": None,
    "observing_strategy.dist_kpc_per_arcsec": 0.3,
    "observing_strategy.inc_deg": 30,
    "observing_strategy.twist_deg": 0,
    "observing_strategy.pointing_deg": "0,0",
    "observing_strategy.pointing_kpc": "",
    "observing_strategy.blur": True,
    "observing_strategy.fwhm": 0.6,
    "observing_strategy.psf": "Gaussian",
    "telescope.type": "SAMI",
    "telescope.fov": None,
    "telescope.aperture_shape": None,
    "telescope.wave_range": "",
    "telescope.wave_centre": None,
    "telescope.wave_res": None,
    "telescope.spatial_res": None,
    "telescope.filter": "r",
    "telescope.lsf_fwhm": None,
    "telescope.signal_to_noise": None
}
```

As before, we can check the output of this request using the code below:

```python
post_build_datacube_request.json()

{'id': 'ee05849c-f89e-4a0b-86a0-66000bfda619',
 'url': 'https://simspin.datacentral.org.au/api/build_datacube/ee05849c-f89e-4a0b-86a0-66000bfda619/',
 'created': '2023-07-26T10:01:10.889497Z',
 'email': '',
 'simspin_file': 'https://simspin.datacentral.org.au/storage/public/build_datacube/ee05849c-f89e-4a0b-86a0-66000bfda619/uploads/SimSpin_example_EAGLE_BC03lr.Rdata',
 'simspin_file_id': None,
 'method': 'velocity',
 'object_name': 'EAGLE_example_SAMI',
 'telescope_name': 'SAMI',
 'observer_name': 'SimSpin',
 'split_save': True,
 'mass_flag': False,
 'observing_strategy': {'dist_z': None,
                        'dist_Mpc': None,
                        'dist_kpc_per_arcsec': 0.3,
                        'inc_deg': 30.0,
                        'twist_deg': 0.0,
                        'pointing_deg': '0,0',
                        'pointing_kpc': '',
                        'blur': True,
                        'fwhm': 0.6,
                        'psf': 'Gaussian'},
 'telescope':  {'type': 'SAMI',
                        'fov': None,
                        'aperture_shape': '',
                        'wave_range': '',
                        'wave_centre': None,
                        'wave_res': None,
                        'spatial_res': None,
                        'filter': 'r',
                        'lsf_fwhm': None,
                        'signal_to_noise': None}}

```

Navigating to the [url link](https://simspin.datacentral.org.au/api/build_datacube/8cec3fe1-571e-4f00-aa6a-e567a057cb60/) above, we can inspect the output of the code within the web interface.

<img align="centre" src="assets/images/build_datacube_api_example_web_interface.png" width="750" height="237" />

Here, we can see that the post has successfully returned a datacube split into a number of FITS files that can be downloaded directly. Alternatively, we can navigate to the web application and view these files interactively. 

To inspect these files and view our output data cube within the web interface, simply change `/api/` in the URL returned by the API to `/app/`.
{: .note} 

<img align="centre" src="assets/images/build_datacube_api_example_app.png" width="750" height="237" />

This observation can be found at the link below. You can modify aspects of the visualisation using the web application options, such as the scale range of the images and the colour mapping. In the application, you can use the "Download Bundle" button on the right to download all of the FITS files produced by the function as a zipped directory.  

<span class="fs-4">
[See the files](https://simspin.datacentral.org.au/app/build_datacube/ee05849c-f89e-4a0b-86a0-66000bfda619/){: .btn .btn-purple }
</span>

---

This concludes the description of using the API directly through Python. You can now produce your own mock observations using SimSpin, without the need to run R locally. These FITS files can then be analysed using your favourite coding language. 

For further examples of using SimSpin or the products produced, check out the rest of our walk-throughs [here](/SimSpin/examples/examples.markdown). 


