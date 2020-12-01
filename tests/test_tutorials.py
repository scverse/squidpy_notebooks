from typing import final

from conftest import MetaTester


class TutorialsMeta(MetaTester):
    pass


@final
class TestTutorials(metaclass=TutorialsMeta):
    pass
