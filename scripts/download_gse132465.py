"""Download GSE132465 UMI count matrix with resume support."""
import urllib.request, ssl, os, time

url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz'
out = r'/path/to/xelox_project\data\scrna\GSE132465_raw_UMI_count_matrix.txt.gz'

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

# Check what we already have
downloaded = os.path.getsize(out) if os.path.exists(out) else 0
print(f"Already downloaded: {downloaded/1024/1024:.1f} MB")

if downloaded > 0:
    # Try to resume
    req = urllib.request.Request(url)
    req.headers['Range'] = f'bytes={downloaded}-'
    try:
        resp = urllib.request.urlopen(req, context=ctx, timeout=30)
        total_expected = int(resp.headers.get('Content-Length', 0)) + downloaded
        print(f"Resuming from {downloaded/1024/1024:.1f} MB, total expected: {total_expected/1024/1024:.1f} MB")
        with open(out, 'ab') as f:
            while True:
                chunk = resp.read(1024*1024)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                print(f"  {downloaded/1024/1024:.1f} MB / {total_expected/1024/1024:.1f} MB ({downloaded/total_expected*100:.0f}%)")
    except Exception as e:
        print(f"Resume failed: {e}")
else:
    # Fresh download
    try:
        req = urllib.request.Request(url)
        resp = urllib.request.urlopen(req, context=ctx, timeout=600)
        total_expected = int(resp.headers.get('Content-Length', 0))
        print(f"New download: {total_expected/1024/1024:.1f} MB")
        with open(out, 'wb') as f:
            while True:
                chunk = resp.read(1024*1024)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                print(f"  {downloaded/1024/1024:.1f} MB / {total_expected/1024/1024:.1f} MB ({downloaded/total_expected*100:.0f}%)")
    except Exception as e:
        print(f"Download failed: {e}")

print(f"Final: {os.path.getsize(out)/1024/1024:.1f} MB" if os.path.exists(out) else "FAILED")
