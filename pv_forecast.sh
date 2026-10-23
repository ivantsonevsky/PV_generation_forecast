#!/bin/bash
# Run the PV forecast and publish the plot through the same GitHub repository
# as the satellite images (see publish_sat.sh).
set -euo pipefail

PY=/Users/ivantsonevsky/opt/anaconda3/bin/python3
WORKDIR=/Users/ivantsonevsky/PV_generation_forecast

SAT_DIR="/Users/ivantsonevsky/SAT/images/latest"       # the folder publish_sat.sh turns into the repo
PUBLISH_SCRIPT="/Users/ivantsonevsky/SAT/publish_sat.sh"   # adjust to where publish_sat.sh lives
PNG="${WORKDIR}/pv_forecast.png"

cd "$WORKDIR"

$PY "${WORKDIR}/pv_forecast_ifs.py" --arrays 3 --capacity 3.5 1.8 3.5 --tilt 30 30 30 \
      --azimuth 65 110 245 --name NE SE SW --inverter-kw 9 --out "$PNG"

# Hand the plot to the publisher: copy under a temporary name, then rename, so
# publish_sat.sh never sees a half-copied file. -p keeps the plot's own
# modification time, which publish_sat.sh records in latest.json.
cp -p "$PNG" "${SAT_DIR}/.pv_forecast.png.tmp"
mv "${SAT_DIR}/.pv_forecast.png.tmp" "${SAT_DIR}/pv_forecast.png"

# publish_sat.sh waits until a file has been unchanged for SETTLE seconds (15 s).
# If it is also run from launchd, you can drop the next block and let that run
# pick the plot up.
if [ -x "$PUBLISH_SCRIPT" ]; then
  sleep 16
  "$PUBLISH_SCRIPT"
else
  echo "Plot copied to ${SAT_DIR}; $PUBLISH_SCRIPT not found, it will be published by the next scheduled run"
fi
