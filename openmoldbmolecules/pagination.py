def get_page_wo_page(currentpath):
    string = currentpath.split('&')
    new = ""
    for item in string:
        if "page" in item:
            pass
        else:
            new += item + "&"
    return new
    