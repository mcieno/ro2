# ro2

Operations Research Homeworks.

University of Padova - Department of Information Engineering (DEI).

## Motivation

Implementing a TSP solver for the "operations research 2" class, a.y. 2019-2020.

## Requirements

To build and run this project you'll need IBM ILOG CPLEX installed on your system.

## Build

```sh
make check  # Ensure library paths are correct
make all    # Build the project
make clean  # Delete build files
```

By default, the project binary will be in `./bin/tsp`.

## Usage

As simple as

```sh
tsp [-?V] [--plot] [--help] [--usage] [--version] TSP_FILE
```

More options are available. Give it a `tsp --help`.

```txt
Usage: tsp [OPTION...] TSP_FILE
Solve a Traveling Salesman Problem instance.

  -c, --cutoff=VALUE         Master cutoff value.
  -j, --threads[=N]          Use multithread. Default 4.
  -m, --memory=SIZE          Available memory (size in MB).
      --name=TSPNAME         Name to assign to this problem.
      --plot                 Draw solution (requires GnuPlot).
  -s, --solver=SOLVER        Solving technique. Default dummy_cplex.
  -t, --timelimit=SECONDS    Maximum time the program may run.
      --verbose, --debug, --trace
                             Set program logging level.
  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version

Mandatory or optional arguments to long options are also mandatory or optional
for any corresponding short options.

Report bugs to {marco.cieno, francesco.cazzaro}@studenti.unipd.it.
```

Sample files are provided in `data`.

Here's an example, just for fun:

```sh
tsp data/att48.tsp --plot --solver=dummy_cplex
```

![Instance plot](assets/att48.instance.png)
![Solution plot](assets/att48.subtour.png)
