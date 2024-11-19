class Primer:
    def __init__(self, sequence, regex, name):
        self.sequence = sequence
        self.regex = regex
        self.name = name
        self.hasForwardMatch = False
        self.hasReverseMatch = False


class RestrictionSite:
    def __init__(self, label, position):
        self.label: str = label
        self.position: int = position
        self._cuts: int = 0

    def set_cut_count(self, count):
        self.cuts = count


class RestrictionSiteCollection:
    def __init__(self):
        self.sites = []

    def add_site(self, site):
        self.sites.append(site)

    def sort(self):
        self.sites.sort(key=lambda x: x.position)
