#!/bin/bash
set -euo pipefail
set +m   # NEW: ensure background jobs stay in same job-control context

# Determine the root directory based on the script location
ROOT=$(dirname "$(readlink -f "$0")")

# Directory to clear before running
CLEAR_DIR="$ROOT/carbonara_runs/smarcalAlpha/fitdata"

# ========= NEW: bookkeeping for clean shutdown =========
PIDS=()
WATCHER_PID=""
cleanup() {
    echo
    echo ">>> Stopping Carbonara frontend..."

    # Stop watcher first
    if [[ -n "${WATCHER_PID}" ]]; then
        kill -INT "$WATCHER_PID" 2>/dev/null || true
        sleep 0.5
        kill -TERM "$WATCHER_PID" 2>/dev/null || true
        sleep 0.5
        kill -KILL "$WATCHER_PID" 2>/dev/null || true
    fi

    # Stop all predictor processes explicitly
    if ((${#PIDS[@]})); then
        echo ">>> Stopping predictor processes (${#PIDS[@]})"
        kill -INT  "${PIDS[@]}" 2>/dev/null || true
        sleep 1
        kill -TERM "${PIDS[@]}" 2>/dev/null || true
        sleep 1
        kill -KILL "${PIDS[@]}" 2>/dev/null || true
    fi

    echo ">>> Stopped."
    exit 0
}
trap cleanup SIGINT SIGTERM
# ======================================================

# Clear the directory
echo "Clearing directory: $CLEAR_DIR"
rm -r "$CLEAR_DIR"/*
mkdir -p "$CLEAR_DIR"

### argv[ 1] scattering data file
ScatterFile=$ROOT/carbonara_runs/smarcalAlpha/Saxs.dat
### argv[ 2] sequence file location
fileLocs=$ROOT/carbonara_runs/smarcalAlpha/
### argv[ 3] restart tag (use to start from existing prediction)
initialCoordsFile=frompdb
### argv[ 4] paired distances file (can be empty)
pairedPredictions=False
### argv[ 5] fixed sections file (again can be empty)
fixedsections=$ROOT/carbonara_runs/smarcalAlpha/varyingSectionSecondary1.dat
### argv[ 6] number of structures
noStructures=1
### argv[ 7] request to apply hydrophobic covering WITHIN monomers will be a list of sections on which to apply it -- Currently not used
withinMonomerHydroCover=none
### argv[ 8] kmin
kmin=0.01
### argv[ 9] kmax
kmax=0.2
### argv[ 10] kmax Start
kmaxStart=0.2
### argv[11] Max number of fitting steps
maxNoFitSteps=10000
### argv[12] prediction file - mol[i] in the fitting folder
predictionFile=$ROOT/carbonara_runs/smarcalAlpha/fitdata
### argv[13] scattering output file
scatterOut=$ROOT/carbonara_runs/smarcalAlpha/fitdata
### argv[14] mixture list file, alist of sets of numbers indicatig the allowed set of mixture percentages of each species (e.g. dimer 20 monomer 80)
mixtureFile=$ROOT/carbonara_runs/smarcalAlpha/mixtureFile.dat
### argv[15] previous fit string in form fitname/mol6Substep_10_1.dat+fitname/mol6Substep_10_2.dat
prevFitStr=$ROOT/carbonara_runs/smarcalAlpha/redundant
### argv[16] log file location
logLoc=$ROOT/carbonara_runs/smarcalAlpha/fitdata
### argv[17] last line of the previous fit log, this is only used for a restart if argv[3] = True
endLinePrevLog=null
### argv[18] is true if we want to apply affine rotations,false if not.
affineTrans=False
### argv[19] is true if we want to use errors in the scattering calculation false if not.
useErrors=True   # IMPORTANT: must be True

# ========= NEW: start watcher (background) =========
WATCHER_SCRIPT="$ROOT/watch_and_backmap.py"
BACKMAP_SCRIPT="$ROOT/backmap_cli.py"
WATCHER_LOG="$predictionFile/watcher.out"

python "$WATCHER_SCRIPT" \
    --watch-dir "$predictionFile" \
    --scenario-root "$ROOT/carbonara_runs/smarcalAlpha" \
    --backmap-script "$BACKMAP_SCRIPT" \
    --max-backmap 3 \
    --defer-backmap-seconds 600 \
    --do-foxs \
    --foxs-py "$ROOT/pyFoXS/pyFoXS/foxs.py" \
    --saxs "$ScatterFile" \
    --max-q "$kmaxStart" \
    > "$WATCHER_LOG" 2>&1 &
WATCHER_PID=$!
echo "Watcher started (PID=$WATCHER_PID)"
# ==================================================

for i in {1..20}
do
    echo "\n"
    echo " >> Run number : $i "
    echo "\n"
    echo "Max number of fitting steps: " $maxNoFitSteps
    echo "\n"

    # NEW: stdbuf ensures live output, PID captured
    stdbuf -oL -eL \
    $ROOT/build/bin/predictStructureQvary \
        $ScatterFile \
        $fileLocs \
        $initialCoordsFile \
        $pairedPredictions \
        $fixedsections \
        $noStructures \
        $withinMonomerHydroCover \
        $kmin \
        $kmax \
        $kmaxStart \
        $maxNoFitSteps \
        $predictionFile/mol$i \
        $scatterOut/scatter$i.dat \
        $mixtureFile \
        $prevFitStr \
        $logLoc/fitLog$i.dat \
        $endLinePrevLog \
        $affineTrans \
        $useErrors \
        > "$predictionFile/run$i.out" 2> "$predictionFile/run$i.err" &

    PIDS+=($!)   # NEW: track PID
done

echo
echo ">>> All runs launched"
echo ">>> Press Ctrl+C to stop everything"
echo

wait
