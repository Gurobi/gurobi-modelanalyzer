from collections import defaultdict
import zipfile
import hashlib
import tempfile
import shutil
import os
import gzip
import bz2


def read_comments(filename):
    # Read up to 50 comment lines from the top of the model file
    # Heuristic: Read all lines on top of the file that start with a non-alphanumeric character
    # Note: OPB and MPS use '*', LP uses '\' for comments

    comment_lines = []
    max_lines = 50
    current_line = 0

    with open(filename, "r") as modelfile:
        for line in modelfile:

            current_line += 1
            line = line.rstrip('\n')

            if len(line) > 0 and not line[0].isalnum():
                comment_lines.append(line)
            else:
                break

            if current_line > max_lines:
                break

    return comment_lines


def process_mps_details(modelfile, data):
    print("Processing MPS file details...")

    blocks = {"NAME", "ROWS", "COLUMNS", "RHS", "RANGES", "BOUNDS", "SOS",
              "QMATRIX", "QSECTION", "QUADOBJ", "QCMATRIX", "ENDATA"}

    blocks_ignored = {"OBJSENSE"}

    line_counter = 0;
    result = defaultdict(list)

    with open(modelfile) as mpsfile:
        for line in mpsfile:
            line_counter += 1
            if line[0] == ' ': continue
            statement = line[:10].upper().rstrip("\n").split(" ")[0].strip()
            if len(statement) > 0:
                result[statement].append(line_counter);

    sections = {}
    other_sections = []

    for section in blocks:
        if result.get(section, None) is None:
            sections[section] = False
        else:
            sections[section] = result[section]

    for section in result.keys():
        if sections.get(section, None) is None and not section in blocks_ignored:
            other_sections.append(section)

    data["fileLines"] = line_counter
    data["mpsSections"] = sections;
    data["mpsOtherSections"] = other_sections;


def process_lp_details(modelfile, data):
    print("Processing LP file details...")

    blocks = [["Minimize", "Minimum", "Min"],
              ["Maximize", "Maximum", "Max"],
              ["Subject To", "Such That", "st.", "st", "s.t.", ],
              ["Lazy Constraints"],
              ["Cones"],
              ["Bounds", "Bound"],
              ["Binaries", "Binary", "Bin"],
              ["Generals", "General", "Gen", "Integers"],
              ["Semi-Continuous", "Semis", "Semi"],
              ["SOS"],
              ["PWLObj"],
              ["End"]];

    max_length = max([len(keyword) for block in blocks for keyword in block]) + 1
    line_counter = 0;

    result = [False for i in range(len(blocks))]

    with open(modelfile) as lpfile:
        for line in lpfile:
            line_counter += 1
            if line[0] == ' ': continue
            statement = line[:max_length];
            for index in range(len(blocks)):
                for keyword in blocks[index]:
                    if statement.lower().startswith(keyword.lower()) and (
                            statement[len(keyword)] == ' ' or statement[len(keyword)] == '\n'):
                        result[index] = [statement[:len(keyword)], line_counter];
                        break

    block_lines = [];

    for i in range(len(blocks)):
        if result[i]:
            block_lines.append(result[i])
        else:
            block_lines.append([blocks[i][0], False])

    data["fileLines"] = line_counter
    data["lpSections"] = block_lines;


def process_opb_details(modelfile, data):
    print("Processing OPB file details...")

    line_counter = 0;

    with open(modelfile) as opbfile:
        for line in opbfile:
            line_counter += 1

    data["fileLines"] = line_counter


def get_file_type(filename):
    mod_exts = {".mps": "MPS", ".rew": "MPS", ".lp": "LP", ".rlp": "LP", ".ilp": "LP", ".opb": "OPB"}
    zip_exts = [".zip", ".gz", ".bz2", ".7z"]

    compression = None
    file_type = None

    filename_lc = filename.lower()

    # Check compression extension and remove from filename if found
    for zip_ext in zip_exts:
        if filename_lc.endswith(zip_ext):
            compression = zip_ext[1:]
            filename_lc = filename_lc[:-len(zip_ext)]
            break

    # Check known extensions
    for mod_ext, ext_ftype in mod_exts.items():
        if filename_lc.endswith(mod_ext):
            file_type = ext_ftype
            break

    print("File Type: %s" % file_type)
    print("Compression: %s" % compression)
    return file_type, compression


def process_uncompressed_model(modelfile, data):
    data["comments"] = read_comments(modelfile)

    # Calculate hash values of uncompressed model
    BLOCKSIZE = 65536
    hasher_md5 = hashlib.md5()
    hasher_sha1 = hashlib.sha1()

    print("Calculating hash values...")
    with open(modelfile, 'rb') as file:
        buf = file.read(BLOCKSIZE)
        while len(buf) > 0:
            hasher_md5.update(buf)
            hasher_sha1.update(buf)
            buf = file.read(BLOCKSIZE)

    data["fileSHA1"] = hasher_sha1.hexdigest()
    data["fileMD5"] = hasher_md5.hexdigest()

    # Read details per model type
    if data["fileType"] == "MPS":
        process_mps_details(modelfile, data)
    elif data["fileType"] == "LP":
        process_lp_details(modelfile, data)
    elif data["fileType"] == "OPB":
        process_opb_details(modelfile, data)


def process_compressed_model(compression, modelfile, data):
    tmpmodelfile_fd, tmpmodelfile_path = tempfile.mkstemp()
    tmpmodelfile = os.fdopen(tmpmodelfile_fd, 'wb')
    tmpmodelfile_archive = None

    try:

        compressed_file = None

        # GZIP, BZIP2, ZIP
        if compression == "gz":
            print("Uncompressing GZIP file...")
            compressed_file = gzip.open(modelfile, 'rb')
        elif compression == "bz2":
            print("Uncompressing BZIP2 file...")
            compressed_file = bz2.BZ2File(modelfile, 'rb')
        elif compression == "zip":
            print("Uncompressing ZIP file...")
            tmpmodelfile_archive = zipfile.ZipFile(modelfile, 'r')
            firstFile = tmpmodelfile_archive.infolist()[0]
            compressed_file = tmpmodelfile_archive.open(firstFile)

        if compressed_file is not None:
            with compressed_file as f_in:
                shutil.copyfileobj(f_in, tmpmodelfile)
                tmpmodelfile.close()
                process_uncompressed_model(tmpmodelfile_path, data)
                return

    finally:

        data["fileSizeUncompressed"] = os.stat(tmpmodelfile_path).st_size

        os.remove(tmpmodelfile_path)
        del tmpmodelfile_archive