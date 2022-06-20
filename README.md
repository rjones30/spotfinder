# spotfinder
**web application for visualization of coherent bremsstrahlung beam properties**

One commonly used method for producing beams of high-energy photons (gamma rays) is to accelerate electrons to high energy and pass them through a material
target. The radiation emitted by the electrons, known as bremsstrahlung, is normally either circularly polarized (if the electrons are longitudinally
polarized) or unpolarized. A third option, linear polarization, is also available if the target ("radiator") is an oriented crystal. The radiation of a
photon from a high-energy electron passing through a crystal is known as "coherent bremsstrahlung", which has been extensively studied by theoretically
and experimentally [1-5]. At sufficiently high electron energies, a coherent bremsstrahlung beam can be produced with nearly ideal properties: high
intensity, monochromatic, and higly polarized. Operating a high-energy photon source that achieves all three of these ideal properties presents a number
of challenges, particularly related to the type and quality of the crystal. The crystal must be as thin as possbile to reduce the effects of multiple
scattering of the electron beam, but thick enough to produce sufficient photon intensity given the limits of the electron beam current imposed by beam
heating and radiation damage to the crystal. The crystal must be a nearly perfect single crystal with good thermal conductivity made from a material
with low atomic number to reduce the background from the incoherent bremsstrahlung process that is always present. Incoherent bremsstrahlung produces a 
low-energy background that is neither monochromatic nor polarized, and degrades the properties of the coherernt bremsstrahlung beam. This combination of
proprties makes **single-crystal diamond** a unique choice for a high-energy coherent bremsstrahlung source. Below is an image of a single-crystal diamond
radiator mounted for use in the coherent bremsstrahlung beamline in Hall D at Jefferson Lab. 
<p align="center">
  <img src="JD70-103_front_view.svg" width="300">
</p>

### X-ray topograhs ###
The above false-color image is an example of what is known as a crystal X-ray topograph. Displayed in the composite graphic is a 7mm x 7mm single-crystal
diamond that has been polished down to a final thickness of just 50 microns (um) and attached n the lower right corner with conductive epoxy to an aluminum
mounting tab. In spite of being so thin, the rigidity and thermal conductivity of diamond is sufficient to make such a target self-supporting and self-cooled
when inserted into a high-energy electron beam. The color in the image represents an elevation map of one of the planes in the diamond crystal, ranging from
low (blue) to high (yellow) relative to an ideal plane aligned with the crystal at three of its corners. The total displacement from low to high in this
image is only 600 nm, but this slight imperfection is enough to significantly distort the properties of the coherent bremsstrahlung beam produced by the
diamond by 12 GeV electrons.

The above image was produced from the combined analysis of two other images known as "X-ray rocking-curve" topographs.
Rocking curves are measured by placing a crystal in an monochromatic X-ray beam and rotating (rocking) the crystal through a very small angular range
around one of the Bragg diffraction peaks of the crystal. As you rotate the crystal around the axes shown in the figure, either &theta;<sub>H</sub> or
&theta;<sub>V</sub>, different parts of the crystal reach a diffraction maximum at slightly different values of the rocking angle because of distortions
of the crystal structure arising from local defects and distributed strain (plastic deformation). Measuring the angle &theta;<sub>H</sub> or
&theta;<sub>V</sub> of the diffraction maximum for every pixel taken with a high-resolution X-ray camera and plotting the result as a false-color image
produces a rocking-curve topograph for one set of crystal planes. A pair of rocking curve topographs, one for &theta;<sub>H</sub> and the other for
&theta;<sub>V</sub> provides complete information about the "tilt" of the crystal planes away from the ideal. One such pair of rocking curves was used
to generate the elevation map of the crystal shown above.

### the beam spot ###
The region of the crystal that is impacted by the high-energy electron beam is herein referred to as *the beam spot*. If the beam spot were point-like
the effects of crystal imperfections would be negligible because a high-quality single-crystal diamond is locally very close to ideal. The effects of
crystal imperfection on the high-energy photon beam come about because the electron beam radiates simultaneously from an extended region of the crystal,
whereas the orientation of the crystal can only be chosen to match the desired kinematics (monochromatic peak energy, polarization) at one point. In this
tool, the electron beam spot is modeled as a two-dimensional gaussian distribution described by an ellipse that defines the one-sigma contour. Three
parameters are needed to define the beam spot ellipse, eg. the major and minor axes and the inclination angle of the major axis. This tool uses a different
but equivalent parameterization: the sigmas of the x and y projections, and the x,y correlation coefficient. Here is a screenshot of the spotfinder tool
invoked with default parameters, showing the beam spot superimposed on the diamnod crystal topograph on the left, and the resulting distribution of tilt
angles &theta;<sub>V</sub> vs. &theta;<sub>H</sub> on the right. Slider controls accompanied by form numeric inputs allow the beam spot shape and position
on the crystal to be varied, while the results in terms of the tilt intensity distribution, the collimated photon beam intensity spectrum, the coherent
enhancement, and the beam linear polarization spectrum are updated in real time.

![image](https://user-images.githubusercontent.com/7832920/174473510-52d7ec20-959a-4d54-aef3-f1d40c54031e.png)
![image](https://user-images.githubusercontent.com/7832920/174476484-a17446f8-7706-4a75-b9e0-b6600ac05a61.png)
![image](https://user-images.githubusercontent.com/7832920/174476497-023c0799-5d79-458b-ae62-87684895a8ff.png)
![image](https://user-images.githubusercontent.com/7832920/174476541-8f9b79ac-c47c-4989-a7d5-c90bf01793d1.png)

## frontend installation ##

Spotfinder is designed to run as a web application served by an apache webserver under mod_wsgi. Installation and setup of an apache webserver is beyond
the scope of this document. In what follows, it is assumed that the webserver is installed and running. First, you need to chose where under your 
DocumentRoot to install spotfinder. In these examples, it is assumed to be /var/www/html/spotfinder, but any other location under DocumentRoot would be
just as good.

```
$ cd /var/www/html
$ sudo git clone https://github.com/rjones30/spotfinder
$ chown -R apache:apache spotfinder
$ cd spotfinder
```

The core component of the spotfinder webapp is the spotfinder.py wsgi script found at the top level in the spotfinder directory you just created.
Open your local copy in a text editor and check that the following lines near the top of the file agree with your local configuration.

```
docroot = "/var/www/html"
topdir = "/spotfinder"
tmpdir = "/spotfinder/tmp"
self_script = "https://my-apache-server/spotfinder"
```

Local scripts running on the webserver must be able to find the spotfinder directory you just created at f"{docroot}{topdir}" and the tmp area where output
graphics files will be stored at f"{docroot}{tmpdir}". Remote users must be able to access files in these same two directories at f"{self_script}{topdir}"
and f"{self_script}{tmpdir}" respectively. Just before the section quoted above, a set of import statements is listed. Review this list and install any
dependent python libraries that are not already installed on your system. You can do this either as the apache user or as root. Most of these are isntalled
with python by default, but some like urllib may need to be downloaded and installed. Do this before continuing.

```
$ sudo pip install shutils urllib2 numpy pickle base64 random
```

Spotfinder depends on the ROOT data analysis and visualization framework from CERN. Pre-built binary packages of ROOT are available at https://root.cern.ch 
for a number of common platforms. ROOT is open-source, so it can also be downloaded and installed from sources using common open-source tools like cmake and
the gnu C++ compilers. If you decide to install ROOT from sources, make sure you enable the python module (pyroot) in your build. Once this is done, update
the spotfinder/setup.sh script to include ROOT in the shell environment. Here is an example that I use on one of my local spotfinder servers. You will need
to replace my local ROOT install path /cvmfs/oasis.opensciencegrid.org/gluex/root-6.22.06/x86_64. For the version of ROOT, any recent release
of ROOT version 6 will work with spotfinder.

```
export PATH=$PATH:/cvmfs/oasis.opensciencegrid.org/gluex/root-6.22.06/x86_64/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/oasis.opensciencegrid.org/gluex/root-6.22.06/x86_64/lib
export PYTHONPATH=$PYTHONPATH:/cvmfs/oasis.opensciencegrid.org/gluex/root-6.22.06/x86_64/lib
```

With ROOT installed and the setup.sh script updated with the ROOT location, you are now ready to complete the installation. As the apache user within the
spotfinder directory you installed under DocumentRoot in a previous step, do the following.

```
$ source setup.sh
$ root -l
root [0] .x makerootvisuals.C
Building the C++ rootvisuals_C.so and supporting libraries
for the spotfinder tool. Make sure that LD_PRELOAD is not
set in your environment, or you are in for a bad surprise!
root [1] .q
$ python3
>>> import spotfinder
>>> ^d
```

If the above steps complete successfully, the frontend installation is nearly finished. All that remains is to configure spotfinder.py as a recognized
wsgi script within your apache server. There is no one right way to do this, but one way that works on a linux system would be to add the following lines
to /etc/httpd/conf.d/ssl.conf:

```
# Setup for osgprod WSGI module [Richard Jones, 11-16-2020]
WSGIScriptAlias /spotfinder /var/www/wsgi-scripts/spotfinder.py
<Directory "/var/www/wsgi-scripts">
# Comment out the following line to enable running of wsgi scripts
# in multiple python subinterpreters, incompatible with matplotlib!
#WSGIApplicationGroup %{GLOBAL}
SetHandler wsgi-script
Options ExecCGI
<IfVersion < 2.4>
   Order allow,deny
   Allow fromall
</IfVersion>
<IfVersion >= 2.4>
   Require all granted
</IfVersion>
</Directory>
```

If you are following closely to the above example, copy the spotfinder.py script that you modified earlier to the new location at /var/www/wsgi-scripts
and assign the ownership to the appropriate user, either apache or root, then restart the apache server, eg. systemctl restart httpd. Direct a web browser
to https://your-apache-server/spotfinder and you should see an image similar to one of the ones above in this README. This completes the frontend installation.

## references ##
[1] G. Diambrini-Palazzi, "High-Energy Bremsstrahlung and Electron Pair Production in Thin Crystals", Revs. Mod. Phys. vol 40 (1968) p. 611.  
[2] U. Timm, "Coherent bremsstrahlung of electrons in crystals", Fortschr. Phys. vol. 17 (1969) p. 765.  
[3] "Coherent Radiation Sources", edited by A. W. Såenz and H. Uberall, Springer-Verlag, Berlin, 1985.  
[4] W. Kaune, G. Miller, W. Oliver, R.W. Williams, and K.K. Young, "Inclusive cross sections for pion and proton production by photons using collimated coherent bremsstrahlung", Phys. Rev. vol 11-3 (1975) pp. 478-494.  
[5] H. Bilokon, G. Bologna, F. Celani, B. D'Ettorre Piazzoli, R. Falcioni, G. Mannocchi, and P. Picchi, "Coherent bremsstrahlung in crystals as a tool for producing high energy photon beams to be used in photoproduction experiments at CERN SPS", Nuclear Inst. and Meth. 204 (1983), pp. 299--310.  
[6] G. Yang, R.T. Jones, F. Klein, K. Finkelstein, K. Livingston, “Rocking Curve Imaging for Diamond Radiator Crystal Selection”, Journal of Diamond & Related Materials 19 (2010) 719.  
[7] K. Finkelstein, R.T. Jones, A. Pauling, D.C. Sagan, Z. Brown, and D S. Misra, “High Resolution, Monochromatic X-ray Topography Capability at CHESS”, Proceedings of the 12th International Conference on Synchrotron Radiation Instrumentation (SRI-2015),  AIP Conf. Proc. 1741(2016) , 010001.  
