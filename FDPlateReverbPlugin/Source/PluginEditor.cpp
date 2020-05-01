/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
FdplateReverbPluginAudioProcessorEditor::FdplateReverbPluginAudioProcessorEditor (FdplateReverbPluginAudioProcessor& p)
    : AudioProcessorEditor (&p), processor (p)
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
    setSize (400, 300);    
 
    // these define the parameters of our slider object
    gainSlider.setSliderStyle (Slider::LinearBarVertical);
    gainSlider.setRange(0.0, 127.0, 1.0);
    gainSlider.setTextBoxStyle (Slider::NoTextBox, false, 90, 0);
    gainSlider.setPopupDisplayEnabled (false, false, this);
    gainSlider.setTextValueSuffix (" Volume");
    gainSlider.setValue(1.0);
 
    // this function adds the slider to the editor
    addAndMakeVisible (&gainSlider);
     
    gainSlideAttach = std::make_unique<AudioProcessorValueTreeState::SliderAttachment>(processor.parameters,
                                                                                       "gain",
                                                                                       gainSlider);
}

FdplateReverbPluginAudioProcessorEditor::~FdplateReverbPluginAudioProcessorEditor()
{
}

//==============================================================================
void FdplateReverbPluginAudioProcessorEditor::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));

    g.setColour (Colours::white);
    g.setFont (15.0f);
    
    g.drawFittedText ("Plate Reverb!", getLocalBounds(), Justification::centred, 1);
}

void FdplateReverbPluginAudioProcessorEditor::resized()
{
    gainSlider.setBounds(40, 30, 30, getHeight() - 60);
}


void FdplateReverbPluginAudioProcessorEditor::sliderValueChanged(Slider* slider)
{
    
}
