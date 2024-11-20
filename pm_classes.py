class Primer:
    def __init__(self, sequence, regex, name) -> None:
        self.sequence: str = sequence
        self.regex: str = regex
        self.name: str = name
        self.hasForwardMatch: bool = False
        self.hasReverseMatch: bool = False


class Match:
    def __init__(self, name, sequence, position) -> None:
        self.name: str = name
        self.sequence: str = sequence
        self.position: int = position
        self.end: int = position + len(sequence)


class MatchCollection:
    def __init__(self) -> None:
        self.matches: list[Match] = []

    def add_match(self, match: Match) -> None:
        self.matches.append(match)

    # TODO might be backwards
    def sort_matches(self) -> None:
        self.matches.sort(key=lambda x: x.position)

    def get_matches(self, start: int, stop: int) -> list[Match]:
        return [match for match in self.matches if start <= match.position <= stop]


class RestrictionSite:
    def __init__(self, label, position) -> None:
        self.label: str = label
        self.position: int = position
        self._cuts: int = 0

    def set_cut_count(self, count: int) -> None:
        self._cuts = count

    def add_cuts(self, count: int) -> None:
        self._cuts += count

    def get_cut_count(self) -> int:
        return self._cuts


class RestrictionSiteCollection:
    def __init__(self) -> None:
        self.sites = []

    def add_site(self, site: RestrictionSite) -> None:
        self.sites.append(site)

    # TODO might be backwards
    def sort_sites(self):
        self.sites.sort(key=lambda x: x.position)
