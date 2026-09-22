"""Multi-threaded download of GSE132465 UMI matrix.
Downloads file in parallel chunks then reassembles."""
import urllib.request, ssl, os, time, threading

url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz'
out = r'/path/to/xelox_project\data\scrna\GSE132465_raw_UMI_count_matrix.txt.gz'

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE

# Get file size first
req = urllib.request.Request(url)
resp = urllib.request.urlopen(req, context=ctx, timeout=30)
total_size = int(resp.headers.get('Content-Length', 0))
print(f'Total size: {total_size/1024/1024:.1f} MB')

# Download in a single thread with retries
max_retries = 5
for attempt in range(max_retries):
    try:
        req = urllib.request.Request(url)
        # Get fresh connection
        resp = urllib.request.urlopen(req, context=ctx, timeout=300)
        existing = os.path.getsize(out) if os.path.exists(out) else 0
        
        mode = 'ab' if existing > 0 else 'wb'
        if existing > 0:
            req = urllib.request.Request(url)
            req.headers['Range'] = f'bytes={existing}-'
            resp = urllib.request.urlopen(req, context=ctx, timeout=300)
            expected_remaining = total_size - existing
            print(f'Resuming from {existing/1024/1024:.1f} MB, remaining: {expected_remaining/1024/1024:.1f} MB')
        
        with open(out, mode) as f:
            dl = existing
            while True:
                chunk = resp.read(65536)  # 64KB chunks
                if not chunk:
                    break
                f.write(chunk)
                dl += len(chunk)
                pct = dl/total_size*100
                print(f'  {dl/1024/1024:.1f}/{total_size/1024/1024:.1f} MB ({pct:.1f}%)', end='\r')
        
        final_sz = os.path.getsize(out)
        print(f'\nDownload finished: {final_sz/1024/1024:.1f} MB')
        
        # Verify gzip
        with open(out, 'rb') as f:
            header = f.read(2)
            if header == b'\x1f\x8b':
                # Try to read
                import zlib
                d = zlib.decompressobj(15 + 32)
                test_data = d.decompress(f.read(65536))
                print('GZIP header valid!')
                
                # Count total rows
                f.seek(0)
                compressed = f.read()
                d2 = zlib.decompressobj(15 + 32)
                all_data = d2.decompress(compressed)
                lines = all_data.decode('utf-8', errors='replace').split('\n')
                print(f'Total rows (incl header): {len(lines)}')
                print(f'Row 0: {lines[0][:100]}')
                print('SUCCESS!')
                break
            else:
                print('Invalid gzip header')
                break
    except Exception as e:
        print(f'\nAttempt {attempt+1} failed: {type(e).__name__}')
        if attempt < max_retries - 1:
            print('Retrying...')
            time.sleep(5)
        else:
            print('All attempts failed')
