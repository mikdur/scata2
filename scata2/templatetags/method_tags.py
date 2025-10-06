from django import template
from django.template.loader import get_template

register = template.Library()


@register.simple_tag(takes_context=True)
def method_results(context):
    method_name = context.get("object").method
    template = get_template("methods/{m}/result.html".format(m=method_name))
    return template.render(context=context.flatten())
