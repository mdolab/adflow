from importlib import import_module
from docutils.parsers.rst.directives.misc import Include
import yaml
import os

# This is in the file scope rather than the class so it can be imported in the conf.py file
TEMP_FILE = "tmp_cost.rst"


# Helper classes to mock the solver class and attributes
class Empty:
    def __del__(self):
        return


class ReturnsOne:
    def __getattribute__(self, name):
        return 1


class CostFunctionsList(Include):
    """
    This is a variant of the OptionsList class in sphinx_mdolab_theme.
    It is an include directive that creates a list of all cost functions and their descriptions.
    """

    optional_arguments = 1
    option_spec = {
        "filename": str,
    }

    # default file name
    filename = "costFunctions.yaml"
    # const strings
    INDENT = "   "
    # class attributes
    module_path = None
    member_name = None
    costFunctions = None
    yaml = None

    def run(self):
        # Set the self.costFunctions attribute
        self.get_costFunctions()
        # Read the descriptions from YAML
        self.get_descriptions()
        # Generate the temporary RST file
        self.generate_temp_file()
        # Reset the self.arguments to just the temp file name for the inherited include directive
        self.arguments = [TEMP_FILE]
        # Actually run the include directive
        return super().run()

    def generate_temp_file(self):
        with open(TEMP_FILE, "w") as f:
            for key in self.costFunctions.keys():
                # The key is the name of the function
                # Write this as data
                f.writelines(f".. data:: {key}")
                f.write("\n")

                # Get the description from the yaml file
                try:
                    desc = self.yaml[key]["desc"]
                except KeyError as e:
                    raise KeyError(
                        f"The description for function '{key}' is missing from the YAML file {self.filename}"
                    ) from e

                # Split the description because it needs to be indented in the RST file
                if "\n" in desc:
                    desc = desc.splitlines()
                else:
                    desc = [desc]

                # Write the description
                f.write(f"\n{self.INDENT}")
                f.writelines(f"\n{self.INDENT}".join(desc))

                # Insert newline before the next data directive
                f.write("\n\n")

    def get_costFunctions(self):
        # Access the class name
        self.module_path, self.member_name = self.arguments[0].rsplit(".", 1)

        # Import the class
        cls = getattr(import_module(self.module_path), self.member_name)

        # Create a mock instance of the class
        inst = Empty()
        inst.adflow = Empty()
        inst.adflow.adjointvars = ReturnsOne()
        inst.adflow.constants = ReturnsOne()

        # Call the private function with the mock class to get the cost functions
        self.costFunctions = cls._getObjectivesAndDVs(inst)[2]

    def get_descriptions(self):
        if "filename" in self.options:
            self.filename = self.options["filename"]
        source_file = self.state.document.attributes["source"]
        self.filename = os.path.join(os.path.dirname(source_file), self.filename)
        if not os.path.isfile(self.filename):
            raise FileNotFoundError(f"The file {self.filename} must exist! Failed module is {self.member_name}.")

        with open(self.filename) as f:
            self.yaml = yaml.load(f, Loader=yaml.FullLoader)


def setup(app):
    app.add_directive("costfunctionslist", CostFunctionsList)
