# Convert the paper's Supplementary Tables (xlsx) to csv for cross-checking (xlsx is large, so it is parsed as a stream with the standard library)
import zipfile, csv, re, sys
import xml.etree.ElementTree as ET
ns = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
src, out = sys.argv[1], sys.argv[2]
z = zipfile.ZipFile(src)
ss = [''.join(t.text or '' for t in si.iter(ns + 't')) for si in ET.fromstring(z.read('xl/sharedStrings.xml')).iter(ns + 'si')]
names = re.findall(r'<sheet name="([^"]+)" sheetId="\d+" r:id="rId(\d+)"', z.read('xl/workbook.xml').decode())
def col(ref):
    n = 0
    for ch in re.match(r'[A-Z]+', ref).group():
        n = n * 26 + ord(ch) - 64
    return n - 1
for name, rid in names:
    fn = f"{out}/{name.replace(' ', '_')}.csv"
    with z.open(f'xl/worksheets/sheet{rid}.xml') as f, open(fn, 'w', newline='') as o:
        w = csv.writer(o)
        for ev, el in ET.iterparse(f):
            if el.tag != ns + 'row':
                continue
            row = {}
            for c in el.iter(ns + 'c'):
                v = c.find(ns + 'v')
                if v is not None:
                    row[col(c.get('r'))] = ss[int(v.text)] if c.get('t') == 's' else v.text
            w.writerow([row.get(i, '') for i in range(max(row) + 1)] if row else [])
            el.clear()
    print(name, '->', fn)
