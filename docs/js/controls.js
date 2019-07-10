function setValue(input, element_id)
{
  if (arguments.length === 2)
  {
    document.getElementById(element_id).innerHTML = input.value;
  }
  switch (input.id)
  {
    case 'amp':
      forceAmp = input.value;
      break;
    case 'freq':
      setForce(input.value)
      break;
    case "display_mode":
      display_3d = input.checked;

      break;
    case "input_mode":
      input_sine = input.checked;
      break;
  }
}
