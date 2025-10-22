Welcome to the IMD-streaming workshop at the MDAnalysis UGM 2025! 
This workshop is a part of the broader IMD Streaming Session. The Streaming Session agenda has also been linked below for reference.
This repository contains all the materials, scripts and instructions for the workshop

If you're interested in using our tools after the workshop and watching a recording, see [post workshop](#post-workshop)

##### Table of Contents  

- [Session agenda](#session-agenda)
- [Streaming Workshop](#streaming-workshop)
  - [1. Codespace environment setup](#1-codespace-environment-setup)
    - [i. Github codespaces in the browser (recommended)](#i-github-codespaces-in-the-browser-recommended)
    - [ii. Github codespaces tunnel from your IDE (VSCode and Pycharm)](#ii-github-codespaces-tunnel-from-your-ide-vscode-and-pycharm)
      - [VSCode](#vscode)
      - [Pycharm](#pycharm)
    - [iii. Local codespace in IDE (VSCode only) (slow, not recommend)](#iii-local-codespace-in-ide-vscode-only-slow-not-recommend)
  - [2. Getting started with the workshop](#2-getting-started-with-the-workshop)
- [Post-workshop](#post-workshop)
  - [Workshop recording and materials](#workshop-recording-and-materials)
  - [Using IMDv3-compatible simulation engines](#using-imdv3-compatible-simulation-engines)
    - [GROMACS](#gromacs)
      - [With docker](#with-docker)
      - [From source](#from-source)
      - [New options](#new-options)
    - [LAMMPS](#lammps)
      - [With docker](#with-docker-1)
      - [From source](#from-source-1)
      - [New options](#new-options-1)
    - [NAMD](#namd)


# Session agenda

If you'd like to follow along with the speakers or use the demo codes after the workshop, all workshop demo code is available in this repo. 

| Time | Topic | Location | Speaker | Code | Presentation
| --- | --- | --- | --- | --- | --- |
| 2:30 PM - 2:50 PM | 🖼️ IMD Streaming Introduction | PSF 186 | Matthias Heyden | | [01-Streaming_Big_Picture-Heyden.pdf](presentations/01-Streaming_Big_Picture-Heyden.pdf)
| 2:50 PM - 3:10 PM | 👀📦 IMDv3 Streaming: Theory, Implementation, Technical Details | PSF 186 | Lawson Woods | | [03-Streaming_MDAnalysis_Functionality-Woods.ipynb](presentations/03-Streaming_MDAnalysis_Functionality-Woods.ipynb)
| 3:10 PM - 3:30 PM | 🚀 IMDv3 in Practice: MD Packages, Performance | PSF 186 | Amruthesh Thirumalaiswamy | | [02-Streaming_MD_Packages_and_IMDClient-Thirumalaiswamy.pdf](presentations/02-Streaming_MD_Packages_and_IMDClient-Thirumalaiswamy.pdf)
| 3:30 PM - 3:45 PM | ☕ Tea / Coffee Break | PSF foyer | | |
| 3:45 PM - 4:05 PM | 👀 Streaming Applications Demo (2 Examples) | PSF 186 | Heekun Cho | [gromacs-demos/vdos/demo.ipynb](gromacs-demos/vdos/demo.ipynb), [namd-demos/ion-flux/ion-flux.ipynb](namd-demos/ion-flux/ion-flux.ipynb) | [04-Application_Velocity_correlation_functions_and_2PT-Cho.pdf](presentations/04-Application_Velocity_correlation_functions_and_2PT-Cho.pdf), [04-Application_Ion_channel_permeation-Cho.pdf](presentations/04-Application_Ion_channel_permeation-Cho.pdf)
| 4:05 PM - 5:05 PM | 🎯 Streaming Workshop | PSF 186 | Amruthesh Thirumalaiswamy | [workshop.ipynb](workshop/workshop.ipynb) |
| 5:05 PM - 5:25 PM | 🌊 Integrating MDAnalysis Streaming Analysis with WESTPA | PSF 186 | Jamie Rowe | |
| 5:25 PM - 5:35 PM | ☕ Tea / Coffee Break | PSF foyer | | |
| 5:35 PM - 5:55 PM | 👀 Streaming Applications Demo | PSF 186 | Heekun Cho | [namd-demos/rmsd-rdf/rmsd-rdf.ipynb](namd-demos/rmsd-rdf/rmsd-rdf.ipynb) |
| 5:55 PM - 6:00 PM | 🚪 Closing Remarks | PSF 186 | Irfan Alibay | | [05-Future_Directions-Heyden.pdf](presentations/05-Future_Directions-Heyden.pdf)


# Streaming Workshop

To get started, we recommend using VSCode in the browser with the Github codespace we've provided which includes all the tools you'll need to get started with live simulation streaming.

## 1. Codespace environment setup

### i. Github codespaces in the browser (recommended)

The easiest way is to simply use this repository to create a codespace.
A workshop environment will be created and VSCode will automatically run in your browser.

Duplicate this tab so you will still have access
to these instructions when the codespace is launched.

Select the green "Code" button and then create a codespace:

![alt text](.media/browser_1.png)
![alt text](.media/browser_2.png)

You're done! The codespace will launch in the current tab. Move on to section 2 to get started with the [workshop](#2-getting-started-with-the-workshop).

### ii. Github codespaces tunnel from your IDE (VSCode and Pycharm)

You can use your own IDE to spin up and connect to a codespace (which GitHub will host). 

#### VSCode

If you have VSCode installed, you can install the 
[codespace extension](https://marketplace.visualstudio.com/items?itemName=GitHub.codespaces). 

After installing, you'll see the "remote explorer" icon on the left.
Sign in if you aren't already.

![alt text](.media/ide_1.png)

Select the dropdown arrow to select "Github codespaces" and
then select the "+" to create a new codespace.

![alt text](.media/ide_2.png)

A dialog will appear. For the repository, enter "amruthesht/imd-workshop-2025". For the branch, select "main"
For the machine type, select "2 cores, 8GB RAM, 32 GB storage"

After that, VSCode will automatically launch a new window which is executing in the codespace workshop environment.
To troubleshoot, see the documentation [here](https://docs.github.com/en/codespaces/developing-in-a-codespace/using-github-codespaces-in-visual-studio-code).

#### Pycharm

A codespace extension is also available for [Pycharm](https://plugins.jetbrains.com/plugin/20060-github-codespaces).

### iii. Local codespace in IDE (VSCode only) (slow, not recommend)

You can also run the workshop activity locally if you have the [devcontainers VScode extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)
and [docker](https://docs.docker.com/engine/install/) installed. 

After docker is installed & enabled and your user has been added to the docker group, run:
```bash
git clone https://github.com/amruthesht/imd-workshop-2025.git
code imd-workshop-2025
```
In VSCode, enter CTRL+SHIFT+P and type: "Dev Containers: Open Folder in Container..." and select
the root of the cloned repo as the folder path. A new window will open which is executing 
in the workshop activity codespace.

## 2. Getting started with the workshop

First, open the [workshop/workshop.ipynb](workshop/workshop.ipynb) Jupyter notebook from this repo in your codespace environment.

![alt text](.media/codespace_1.png)

Before running any code, click the "Select kernel" button in the upper right corner of the Jupyter notebook.

![alt text](.media/codespace_2.png)

Select "Python environments" and then the "workshop" environment.

![alt text](.media/codespace_3.png)
![alt text](.media/codespace_4.png)

Now you're ready to start the workshop! Follow the instructions in the notebook to complete the streaming activities.

If you need help, you can refer to the [workshop/solutions.ipynb](workshop/solutions.ipynb) notebook which contains solutions/ expected output for all the try-yourself sections.

# Post-workshop

If you are interested in using our tools, please feel free to reach out for support, bug reports, or for sharing your ideas!

The best way to reach us is on the [MDAnalysis Discord](https://discord.gg/fXTSfDJyxE) in the '#streaming' channel. You can also reach out to us via email (workshops@mdanalysis.org)

## Workshop recording and materials

* A recording of the workshop will be made available after the event.
* All workshop materials are available in _this_ repository. Broader Streaming Session materials are linked from the [Session agenda](#session-agenda).
* Below, we provide instructions for using the simulation engines integrated with IMDv3 capability. You can either use *docker images* (for GROMACS and LAMMPS) or *build your own version* (GROMACS, LAMMPS, NAMD).

## Using IMDv3-compatible simulation engines

For docker usage, ensure [docker](https://docs.docker.com/engine/install/) is installed and the [nvidia container toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/index.html) is installed if using GPU builds.

### GROMACS

#### With docker


First, pull the container:

```bash
# CPU-only build
docker pull ghcr.io/becksteinlab/streaming-md-docker:main-common-cpu

# CUDA build
docker pull ghcr.io/becksteinlab/streaming-md-docker:main-common-gpu
```

To run GROMACS, do:

```bash
# CPU
docker run -v /path/to/input/files:/home/conda:rw -p 8889:8889 \
    ghcr.io/becksteinlab/streaming-md-docker:main-common-cpu bash -c 'gmx <cmd>'

# GPU
docker run -v /path/to/input/files:/home/conda:rw --runtime=nvidia --gpus=all -p 8889:8889 \
    ghcr.io/becksteinlab/streaming-md-docker:main-common-gpu bash -c 'gmx <cmd>'
```

#### From source

The modified codes are available in [this GROMACS fork](https://gitlab.com/ljwoods2/gromacs).

First, clone in the repo:
```
git clone https://gitlab.com/ljwoods2/gromacs.git
git checkout imd-v3
```

For build instructions, see the [GROMACS installation instructions](https://manual.gromacs.org/documentation/current/install-guide/index.html)

#### New options

New MDP file options (subject to change as we work with GROMACS developers):
```
IMD-group               = <group> ; Use 'System' to send the entire system via IMD (inherited from IMDv2)
IMD-version             = <2 | 3> ; Defaults to 2 for backwards compatibility
IMD-nst                 = <nst>   ; Number of integration steps between simulation frames communicated via IMD, defaults to 100
IMD-time                = <yes | no> ; Whether to send time and step information via IMD, defaults to 'no'
IMD-box                 = <yes | no> ; Whether to send box dimension information via IMD, defaults to 'no'
IMD-coords              = <yes | no> ; Whether to send atomic coordinate information via IMD, defaults to 'no'
IMD-vels                = <yes | no> ; Whether to send atomic velcities information via IMD, defaults to 'no'
IMD-forces              = <yes | no> ; Whether to send atomic forces information via IMD, defaults to 'no'
IMD-unwrap              = <yes | no> ; Whether to unwrap molecules to make them appear whole, defaults to 'no'
IMD-energies            = <yes | no> ; Whether to send system energy information via IMD, defaults to 'no'
```
Note that new options will not be used if "IMD-version" is set to 2.

`mdrun` command line options for IMD are inherited from IMDv2, see [gmx-mdrun](https://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html)

### LAMMPS

#### With docker

First, pull the container:

```bash
# CPU-only build
docker pull ghcr.io/becksteinlab/streaming-md-docker:main-common-cpu

# CUDA build
docker pull ghcr.io/becksteinlab/streaming-md-docker:main-common-gpu
```

To run LAMMPS, do:

```bash
# CPU
docker run -v /path/to/input/files:/home/conda:rw -p 8889:8889 \
    ghcr.io/becksteinlab/streaming-md-docker:main-common-cpu bash -c 'lmp < </path/to/infile>'

# GPU
docker run -v /path/to/input/files:/home/conda:rw --runtime=nvidia --gpus=all -p 8889:8889 \
    ghcr.io/becksteinlab/streaming-md-docker:main-common-gpu bash -c 'lmp < </path/to/infile>'
```

#### From source

The modified codes are available in [this LAMMPS fork](https://github.com/ljwoods2/lammps).

First, clone in the repo:
```
git clone https://github.com/ljwoods2/lammps.git
git checkout imd-v3-integration
```

Build instructions are available in the [LAMMPS installation instructions](https://docs.lammps.org/Install.html)

#### New options

Original options in the IMD fix are available [here](https://docs.lammps.org/fix_imd.html).

With our modifications:

```
fix ID group-ID imd <imd_port> [trate <imd_trate>] [version (2|3)] [unwrap (on|off)] [fscale <imd_fscale>] [time (on|off)] [box (on|off)] [coordinates (on|off)] [velocities (on|off)] [forces (on|off)]
```

'version' will default to 2 for backward compatibility, in which case the new options (time, box, positions, etc) will have no effect.

### NAMD

Due to restrictions on distributing NAMD, we are unable to provide a pre-built docker image. However, we provide a patch for NAMD 3.0 to enable IMDv3 compatibility.

#### IMDv3 patch

One can register for and download the NAMD 3.0 source code from the [NAMD website](https://www.ks.uiuc.edu/Development/Download/download.cgi?UserID=&AccessCode=&ArchiveID=1712).

The IMDv3 patch is available as a `*.diff` file in this [repository](https://github.com/amruthesht/namd-3.0-IMDv3-patch). To apply the patch, navigate to the root directory of the NAMD source code and run:

```
  cd /path/to/namd-3.0-source-repo
  patch -p1 < /path/to/namd-3_0-IMDv3.diff
```

Once this is done, the source code will be patched with the new IMDv3 protocol. Detailed compile and build instructions can be found in the `IMDv3-dev.md` file in the patched repository.

#### New options

IMD based options/settings can be set in the NAMD input configuration file.

Previously available options for IMD version 2 in NAMD are available [here](https://www.ks.uiuc.edu/Research/namd/3.0/ug/node49.html).

The following new options are available as a part of the IMDv3 protocol:

```bash
# IMD version -- 2 for VMD and 3 for latest protocol, defaults to 2
IMDversion     3
# IMD session info settings
# IMDsendPositions -- sending positions of entire system
IMDsendPositions        yes
# IMDsendEnergies -- sending energies and bonded, non-bonded and other contributions
IMDsendEnergies     yes
# IMDsendTime -- sending time information (time, dt, step)
IMDsendTime        yes
# IMDsendBoxDimensions -- sending box dimensions (lattice vectors a, b, c)
# If box dimensions are not defined, default unit box is sent
IMDsendBoxDimensions       yes
# IMDsendVelocities -- sending velocities of entire system
IMDsendVelocities       yes
# IMDsendForces -- sending forces on all atoms
IMDsendForces      yes
# IMDwrapPositions -- wrapping positions to box; applicable when IMDsendPositions is yes
IMDwrapPositions       yes
```

When `IMDversion` is set to 2, the new options (`IMDsendTime`, `IMDsendBoxDimensions`, `IMDsendVelocities`, `IMDsendForces`, `IMDwrapPositions`) will have no effect.
