from snakemake.script import snakemake
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
import shutil

import sys


def render_jinja2(template_folder: Path, variables: dict, output: Path):
    template_file = Path(variables["input"]["template"])
    shutil.copy(template_file, template_folder)

    # Handle potential_file only if it's specified
    potential_file = variables["input"].get("potential_file")
    if potential_file:
        potential_file = Path(potential_file)
        shutil.copy(potential_file, template_folder)
        variables["input"]["potential_file"] = Path(
            variables["input"]["potential_file"]
        ).name

    env = Environment(loader=FileSystemLoader(str(template_folder)))
    template = env.get_template(template_file.name)

    with open(output, "w") as f:
        f.write(template.render(**variables))


variables = dict(params=dict(snakemake.params), input=dict(snakemake.input))

output = snakemake.output[0]  # the actual rendered file
template_folder = snakemake.params[
    "template_folder"
]  # path to where the template file is. I guess this could be relative?

# Make the folder for the templates
Path(template_folder).mkdir(parents=True, exist_ok=True)

render_jinja2(template_folder, variables, output)
