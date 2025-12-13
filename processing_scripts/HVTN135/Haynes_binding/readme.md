The lab originally uploaded binding data for HVTN115 and 135 in one input file.
We have a script that simultaneously processes that data, and saves outputs. It is stored under the HVTN115 dir in this repo.

The lab later uploaded binding data just for HVTN135.
The script in this dir processes the new upload, reads in our prior outputs from the initial upload, concats them and saves.

12/12: The lab uploaded new data just for HVTN135 infants at visitno 1 and 12. we processed that in "process_12_11_data.py" and concatenated it onto the 11/25 outputs and reprovisioned the data.