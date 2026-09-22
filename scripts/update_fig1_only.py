#!/usr/bin/env python3
"""
Quick update: Replace Figure 1 in the latest manuscript with redrawn version.
"""
import os
from docx import Document
from docx.shared import Inches

BASE = r'/path/to/xelox_project'
SRC = os.path.join(BASE, 'reports', 'manuscript_draft_v5.8_submission.docx')
OUT = os.path.join(BASE, 'reports', 'manuscript_draft_v5.8_submission.docx')
FIG = os.path.join(BASE, 'results', 'figures_png', 'hd', 'Fig1_pipeline.png')

doc = Document(SRC)
print(f'Loaded: {SRC}')
print(f'Paragraphs: {len(doc.paragraphs)}')

# Find and replace Figure 1 image
fig1_para = None
for i, para in enumerate(doc.paragraphs):
    # Find the paragraph where Fig1 was inserted (after "discovery cohort of 164")
    t = para.text
    if 'Figure 1.' in t and ('pipeline' in t.lower() or 'study design' in t.lower()):
        fig1_para = para
        print(f'Found Figure 1 caption: P{i}')
        
        # The image is in the previous paragraph (inserted before caption)
        prev = doc.paragraphs[i - 1]
        for run in prev.runs:
            if run._element.findall('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}drawing'):
                # Clear existing image runs
                for r in prev.runs:
                    r._element.getparent().remove(r._element)
                print('  Cleared old image')
                break
        
        # Insert new image before the caption
        new_para = doc.add_paragraph()
        run = new_para.add_run()
        run.add_picture(FIG, width=Inches(6.0))
        # Move new_para before caption
        caption_elem = para._element
        caption_elem.addprevious(new_para._element)
        print(f'  Inserted new Figure 1')
        break

if fig1_para is None:
    print('[WARN] Figure 1 caption not found. Trying keyword approach...')
    # Fallback: find by keyword
    for i, para in enumerate(doc.paragraphs):
        if 'discovery cohort of 164' in para.text:
            print(f'Found keyword at P{i}')
            # Insert image before this paragraph
            new_para = doc.add_paragraph()
            run = new_para.add_run()
            run.add_picture(FIG, width=Inches(6.0))
            para._element.addprevious(new_para._element)
            print(f'  Inserted new Figure 1')
            break

doc.save(OUT)
print(f'\nSaved: {OUT}')
print('Done!')
