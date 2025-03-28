from scata2.models import ScataFile, ScataTagSet
import re


def parse_tagset(pk):
    tagset = ScataTagSet.objects.get(pk=pk)
    file = tagset.tagset_file.file.open(mode="rt")

    if not file:
        tagset.is_valid = False
        tagset.validated = True
        tagset.errors = "Internal error. File not found."
        tagset.save()
        return "failed"
    
    if file.size > 100000:
        tagset.is_valid = False
        tagset.validated = True
        tagset.errors = "File to big"
        tagset.save()
        return "failed"

    tag_len = 0
    errors = ""
    success = True
    line_no = 0
    err_cnt = 0
    tags = dict()

    for line in file:
        if err_cnt > 20:
            errors += "Too many errors, bailing out\n"
            tagset.is_valid = False
            tagset.validated = True
            tagset.errors = errors
            tagset.save()
            return "failed"


        line_no += 1

        line = line.rstrip()
        if line == "":
            continue
        
        s = line.split(";")

        if len(s) < 2:
            continue

        if tag_len == 0:
            tag_len = len(s[1])

        if len(s[1]) != tag_len:
            errors += "Tag on line {line} is not {bp} bp, as tags on previous lines\n".format(
                line = line_no, bp=tag_len)
            err_cnt += 1
            continue

        s[1] = s[1].upper()
            
        if not re.match("^[ACTG]+$", s[1]):
            errors += "Tag on line {line} contains invalid bases\n".format(line=line_no)
            err_cnt += 1
            continue

        if not re.match("^[-A-Za-z0-9_]+$",s[0]):
            errors += "Tag name on line {line} contains invalid characters\n".format(line=line_no)
            err_cnt += 1
            continue

        if len(s) > 2:
            for i in range(2, len(s)):
                if not re.match("^[-A-Za-z0-9_]+$", s[i]):
                    errors += "Pairing tag name {num} on line {line} "
                    "contains invalid characters\n".format(name=i + 1, line=line_no)
                    err_cnt += 1
                    continue
                    
        tags[s[1]]={"name": s[0],
                    "mates": set(s[2:])}



    tagset.validated = True
    tagset.is_valid = True
    tagset.errors = errors
    tagset.tags = tags
    tagset.num_tags = len(tags)
    tagset.save()

    return "success"
