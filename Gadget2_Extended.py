from amuse.community.gadget2.interface import Gadget2


class Gadget2_Extended(Gadget2):
    def __init__(self, **options):
        super().__init__(**options)
        self.__init__()

