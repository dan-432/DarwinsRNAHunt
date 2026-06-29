import sys
import argparse
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas

def parse_stockholm(sto_file):
    sequences = []
    structure = None

    with open(sto_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('#=GC SS_cons'):
                structure = line.split(None, 2)[2] if len(line.split(None, 2)) > 2 else ""
            elif line.startswith('#') or line.startswith('//') or not line.strip():
                continue
            elif ' ' in line:
                parts = line.split(None, 1)
                if len(parts) == 2:
                    name, seq = parts
                    sequences.append((name, seq))

    return sequences, structure

def create_visualization_pdf(sto_file, pdf_file, motif_id):
    """Create a simple PDF visualization of the structure."""
    sequences, structure = parse_stockholm(sto_file)

    c = canvas.Canvas(pdf_file, pagesize=letter)
    width, height = letter

    c.setFont("Helvetica-Bold", 16)
    c.drawString(1*inch, height - 1*inch, f"RNA Structure Motif: {motif_id}")

    c.setFont("Helvetica", 10)
    c.drawString(1*inch, height - 1.3*inch, f"Number of sequences: {len(sequences)}")

    if structure:
        c.setFont("Helvetica-Bold", 12)
        c.drawString(1*inch, height - 1.7*inch, "Consensus Secondary Structure:")
        c.setFont("Courier", 8)

        y_pos = height - 2*inch
        chunk_size = 80
        for i in range(0, len(structure), chunk_size):
            chunk = structure[i:i+chunk_size]
            c.drawString(1*inch, y_pos, chunk)
            y_pos -= 0.2*inch

    c.setFont("Helvetica-Bold", 12)
    y_pos = height - 2.5*inch if not structure else y_pos - 0.3*inch
    c.drawString(1*inch, y_pos, "Alignment (first 10 sequences):")

    c.setFont("Courier", 7)
    y_pos -= 0.2*inch

    for name, seq in sequences[:10]:
        if y_pos < 1*inch:
            c.showPage()
            y_pos = height - 1*inch

        display_name = name[:30] + "..." if len(name) > 30 else name
        display_seq = seq[:80] + "..." if len(seq) > 80 else seq

        c.drawString(1*inch, y_pos, f"{display_name}")
        y_pos -= 0.15*inch
        c.drawString(1*inch, y_pos, f"  {display_seq}")
        y_pos -= 0.25*inch

    if len(sequences) > 10:
        c.drawString(1*inch, y_pos, f"... and {len(sequences) - 10} more sequences")

    c.save()

def main():
    parser = argparse.ArgumentParser(description='Create PDF visualization of RNA structure motif from Stockholm alignment')
    parser.add_argument('--stockholm', required=True, help='Input Stockholm alignment file')
    parser.add_argument('--output', required=True, help='Output PDF file')
    parser.add_argument('--motif-id', required=True, help='Motif identifier')

    args = parser.parse_args()

    try:
        create_visualization_pdf(args.stockholm, args.output, args.motif_id)
        print(f"PDF visualization created successfully: {args.output}", file=sys.stderr)
    except Exception as e:
        print(f"Failed to create PDF: {e}", file=sys.stderr)
        try:
            c = canvas.Canvas(args.output, pagesize=letter)
            c.setFont("Helvetica", 12)
            c.drawString(100, 700, f"RNA Structure Motif: {args.motif_id}")
            c.drawString(100, 680, "Visualization unavailable - see .sto file for details")
            c.save()
            print(f"Created placeholder PDF: {args.output}", file=sys.stderr)
        except Exception as e2:
            print(f"Failed to create placeholder PDF: {e2}", file=sys.stderr)
            sys.exit(1)

if __name__ == "__main__":
    main()