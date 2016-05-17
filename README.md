# assignment-2

## Requirements

1. beautifulsoup4
2. goatools
  * fisher
    * numpy

## Directory tree

The following directory tree is assumed by `Assignment2.py`:

```
Assignment2
|
├───assignment2 (this Git repository)
|   │   Assignment2.py
|   |   run.sh
|   │   go-basic.obo
|   │   ...
|   |
|
|
└───DmGOs
    ├───Kc167_P53_NT60_A_GO
    |   |   geneOntology.html
    |   |   ...
    |
    ├───Kc167_P53_NT60_B_GO
    |   ...
```

## Setting up directories on HPC

Make sure you have properly configured your HPC account's Git settings AND added your HPC SSH key to your GitHub account. ([wiki](https://github.com/zhoulab/assignment-2/wiki/Using-GitHub-with-HPC))

1. Log in to HPC:

    ```
    ssh <YOUR_HPC_USERNAME>@hipergator.rc.ufl.edu
    ```

2. Make a directory for the project:

    ```
    mkdir Assignment2
    ```

3. Make a directory for the GO folders (we will need to populate it using `scp` outside of the SSH session): 

    ```
    mkdir Assignment2/DmGOs
    ````

4. Clone this git repository (will create under the `Assignment2` directory:

    ```
    git clone git@github.com:zhoulab/assignment-2.git Assignment2/assignment2
    ```

5. Exit the SSH session (`CTRL+D`) and copy an existing DmGOs folder from your local machine into the `DmGOs` directory.

    ```
    scp -r /path/to/DmGOs/ <YOUR_HPC_USERNAME>@hipergator.rc.ufl.edu:Assignment2/DmGOs/
    ```

## Building

Make sure you are in the `assignment2` directory.

Run `Assignment.py` with `python` and `virtualenv`:

```
./run.sh
```

The repository comes with a virtual environment (`ve/`). If `run.sh` does not work, try the teardown and bootstrap process:

```
./teardown.sh
```

```
./bootstrap.sh
```
