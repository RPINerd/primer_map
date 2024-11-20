"""
    Layout classes for drawing the primer mapping
"""

import re

from pm_utils import right_num


class LayoutComponent:
    def __init__(self):
        self.characters = []
        self.positionLabel = 1

    def write_layout_component(self, start, stop):
        # abstract method
        pass

    def set_characters(self, text):
        if "." in text:
            self.characters = list(text)

    def is_room(self, start, stop):
        range_to_check = "".join(self.characters[start:stop])
        return not any(c.isalnum() for c in range_to_check)


class TranslationComponent(LayoutComponent):
    def write_layout_component(self, start, stop):
        print(right_num(self.positionLabel, "", 8, ""))
        print("".join(self.characters[start:stop]) + "\n")
        self.positionLabel += (stop - start) // 3


class UppercaseTranslationComponent(LayoutComponent):
    def write_layout_component(self, start, stop):
        text_to_write = "".join(self.characters[start:stop]) + "\n"
        if any(c.isalnum() for c in text_to_write):
            print(right_num(self.positionLabel, "", 8, ""))
            self.positionLabel += sum(1 for c in text_to_write if c.isupper())
            print(text_to_write)


class DnaComponent(LayoutComponent):
    def write_layout_component(self, start, stop):
        print(right_num(self.positionLabel, "", 8, ""))
        print("".join(self.characters[start:stop]) + "\n")
        self.positionLabel += stop - start


class RulerComponent(LayoutComponent):
    def write_layout_component(self, start, stop):
        print(right_num(self.positionLabel, "", 8, ""))
        text = "".join(self.characters[start:stop])
        text = re.sub(r"^(\d+)", lambda m: " " * len(m.group(1)), text)
        text = re.sub(r"(\d+)$", lambda m: " " * len(m.group(1)), text)
        print(text + "\n")
        self.positionLabel += stop - start

    def build_ruler(self):
        sequence = "".join(self.characters)
        count = 0
        spacing = "         "
        sequence = re.sub(r"(.{1,10})", lambda m: (spacing if count == 0 else str(count) + spacing)[:10], sequence)
        self.characters = list(sequence)
