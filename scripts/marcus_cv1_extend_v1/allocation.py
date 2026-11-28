import datetime, os, subprocess, sys, re, time

def deadline(text):
    match=re.search(r'\bEndTime=([^ ]+)',text)
    assert match and match[1] not in ('Unknown','None'), 'Live scheduler deadline missing'
    return datetime.datetime.fromisoformat(match[1]).timestamp()

if sys.argv[1]=='selftest':
    assert deadline('JobId=1 EndTime=2026-11-28T12:00:00 X=2')==datetime.datetime(2026,11,28,12).timestamp()
    for s in ('EndTime=Unknown','JobId=1'):
        try:deadline(s)
        except AssertionError:pass
        else:raise AssertionError('Invalid deadline accepted')
    live=subprocess.check_output(['scontrol','show','job',os.environ['SLURM_JOB_ID'],'-o'],text=True)
    assert deadline(live)>time.time()
    print('LIVE_ALLOCATION_DEADLINE_PASS')
else:
    live=subprocess.check_output(['scontrol','show','job',os.environ['SLURM_JOB_ID'],'-o'],text=True)
    end=deadline(live);assert end>time.time()+900
    print(end)
