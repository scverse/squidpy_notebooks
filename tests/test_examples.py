from typing import final

from conftest import MetaTester


class ExamplesMeta(MetaTester):
    pass


@final
class TestExamples(metaclass=ExamplesMeta):
    pass
