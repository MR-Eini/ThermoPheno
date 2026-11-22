"""Extract the documented ThermoPheno weather example from the DWD archive."""
import argparse
import csv
from datetime import date, timedelta
import io
from pathlib import Path
from zipfile import ZipFile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("archive", type=Path, help="DWD daily KL historical ZIP for station 04709")
args = parser.parse_args()

start = date(2005, 1, 1)
end = date(2005, 12, 31)
rows = []
with ZipFile(args.archive) as archive:
    (member,) = [name for name in archive.namelist() if name.startswith("produkt_klima_tag_")]
    with io.TextIOWrapper(archive.open(member), encoding="ascii") as stream:
        for row in csv.DictReader(stream, delimiter=";"):
            fields = {key.strip(): value.strip() for key, value in row.items()}
            day = date.fromisoformat(fields["MESS_DATUM"])
            if start <= day <= end:
                if fields["STATIONS_ID"].zfill(5) != "04709":
                    raise ValueError("The archive is not for DWD station 04709.")
                tmin, tmax = fields["TNK"], fields["TXK"]
                if not -80 < float(tmin) <= float(tmax) < 80:
                    raise ValueError(f"Missing or invalid temperatures for {day}.")
                rows.append((day, tmin, tmax, "04709"))

expected_days = [start + timedelta(days=i) for i in range((end - start).days + 1)]
if [row[0] for row in rows] != expected_days:
    raise ValueError("The archive does not supply every day of 2005 in date order.")

destination = Path(__file__).resolve().parents[1] / "inst/extdata/dwd_simmern_wahlbach_2005.csv"
with destination.open("w", newline="", encoding="utf-8") as stream:
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(("date", "tmin", "tmax", "station"))
    writer.writerows(rows)
