"""Refine only uncertified boxes; preserve the first run as a separate file."""
from pathlib import Path
import json
import time
from certify_outside_one_cycle import I, S, certify_box

ROOT = Path(__file__).parent


def main():
    started = time.time()
    report = json.loads((ROOT/'outside_one_cycle_certificate.json').read_text(encoding='utf-8'))
    leaves = [{**row, 'cells': report['cells']} for row in report['leaves']]
    stack = report['unresolved'][:]
    unresolved = []
    visited = 0
    while stack:
        row = stack.pop()
        bl,bh = row['b']
        rl,rh = row['ratio']
        path = row['path']
        visited += 1
        reason,margin = certify_box(I(bl,bh),I(rl,rh),cells=96)
        if reason:
            leaves.append({**row,'reason':reason,'margin':margin,'cells':96})
        elif len(path) >= 24:
            unresolved.append(row)
        elif bh-bl >= 2*(rh-rl):
            mid = (bl+bh)//2
            stack.extend((dict(path=path+'0',b=[bl,mid],ratio=[rl,rh]),
                          dict(path=path+'1',b=[mid,bh],ratio=[rl,rh])))
        else:
            mid = (rl+rh)//2
            stack.extend((dict(path=path+'0',b=[bl,bh],ratio=[rl,mid]),
                          dict(path=path+'1',b=[bl,bh],ratio=[mid,rh])))
    counts = {}
    for row in leaves:
        key = row['reason']+'_'+str(row['cells'])
        counts[key] = counts.get(key,0)+1
    report.update(status='complete' if not unresolved else 'incomplete',
                  leaves=sorted(leaves,key=lambda x:x['path']),unresolved=unresolved,
                  leaf_count=len(leaves),counts=counts,refinement_visited=visited,
                  refinement_seconds=time.time()-started)
    (ROOT/'outside_one_cycle_certificate_refined.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('leaves','unresolved')},indent=2))
    print(json.dumps(unresolved,indent=2))


if __name__ == '__main__':
    main()
