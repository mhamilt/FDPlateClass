/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
FdplateReverbPluginAudioProcessor::FdplateReverbPluginAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
     : AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  AudioChannelSet::stereo(), true)
                      #endif
                       .withOutput ("Output", AudioChannelSet::stereo(), true)
                     #endif
                       ),
#endif
parameters(*this, nullptr, "ParamTreeExample",
{
    std::make_unique<AudioParameterFloat>("wetdry", "Wetdry", NormalisableRange<float> (0.0f, 0.5f), 0.5f)
})
{
    gainParam = parameters.getRawParameterValue("wetdry");
}

FdplateReverbPluginAudioProcessor::~FdplateReverbPluginAudioProcessor()
{
}

//==============================================================================
const String FdplateReverbPluginAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool FdplateReverbPluginAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool FdplateReverbPluginAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool FdplateReverbPluginAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double FdplateReverbPluginAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int FdplateReverbPluginAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int FdplateReverbPluginAudioProcessor::getCurrentProgram()
{
    return 0;
}

void FdplateReverbPluginAudioProcessor::setCurrentProgram (int index)
{
}

const String FdplateReverbPluginAudioProcessor::getProgramName (int index)
{
    return {};
}

void FdplateReverbPluginAudioProcessor::changeProgramName (int index, const String& newName)
{
}

//==============================================================================
void FdplateReverbPluginAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    plateReverb.setup(sampleRate, FDPlate::PlateParameters());
//    plateReverb.printInfo();
}

void FdplateReverbPluginAudioProcessor::releaseResources()
{
    
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool FdplateReverbPluginAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
    return true;
}
#endif

void FdplateReverbPluginAudioProcessor::processBlock (AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{
    ScopedNoDenormals noDenormals;
//    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();

    const float magicNumber = 400.00;
       
    float  wet = cos(*parameters.getRawParameterValue("wetdry") * MathConstants<float>::pi);
    float  dry = sin(*parameters.getRawParameterValue("wetdry") * MathConstants<float>::pi);
    
    auto* outputData = buffer.getWritePointer (0);
    for (int i = 0; i < buffer.getNumSamples(); ++i)
    {
        outputData[i] = (wet * buffer.getReadPointer (0)[i])
        + (dry * (plateReverb.reverb(buffer.getReadPointer (0)[i]) * magicNumber));
    }
    for (int channel = 1; channel < totalNumOutputChannels; ++channel)
    {
        buffer.copyFrom(channel,
                        0,
                        outputData,
                        buffer.getNumSamples());
    }
    
}

//==============================================================================
bool FdplateReverbPluginAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

AudioProcessorEditor* FdplateReverbPluginAudioProcessor::createEditor()
{
    return new FdplateReverbPluginAudioProcessorEditor (*this);
}

//==============================================================================
void FdplateReverbPluginAudioProcessor::getStateInformation (MemoryBlock& destData)
{
//     MemoryOutputStream (destData, true).writeFloat (*gain);
}

void FdplateReverbPluginAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
//    *gain = MemoryInputStream (data, static_cast<size_t> (sizeInBytes), false).readFloat();
}

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new FdplateReverbPluginAudioProcessor();
}
