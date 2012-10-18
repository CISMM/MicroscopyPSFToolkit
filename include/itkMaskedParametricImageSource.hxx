/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkMaskedParametricImageSource_hxx
#define __itkMaskedParametricImageSource_hxx

#include "itkMaskedParametricImageSource.h"

namespace itk
{

template< class TOutputImage >
MaskedParametricImageSource< TOutputImage >
::MaskedParametricImageSource()
{
  m_DelegateImageSource = NULL;
  m_EnabledArray.clear();
}

template< class TOutputImage >
MaskedParametricImageSource< TOutputImage >
::~MaskedParametricImageSource()
{
  m_DelegateImageSource = NULL;
}

template< typename TOutputImage >
void
MaskedParametricImageSource< TOutputImage >
::SetDelegateImageSource(DelegateImageSourceType * delegateSource)
{
  if ( delegateSource != m_DelegateImageSource )
    {
    m_DelegateImageSource = delegateSource;

    // Enable all parameters from the delegate source
    if ( m_DelegateImageSource )
      {
      m_EnabledArray = EnabledArrayType( m_DelegateImageSource->GetNumberOfParameters() );
      for ( unsigned int i = 0; i < m_DelegateImageSource->GetNumberOfParameters(); ++i )
        {
        m_EnabledArray[i] = true;
        }
      }
    else
      {
      m_EnabledArray = EnabledArrayType( 0 );
      }

    this->Modified();
    }
}

template< typename TOutputImage >
void
MaskedParametricImageSource< TOutputImage >
::SetParameters(const ParametersType & parameters)
{
  if ( m_DelegateImageSource.GetPointer() == NULL )
    {
    itkExceptionMacro( << "Cannot set parameters when DelegateImageSource is not set" );
    }

  // Set only the active parameters in the delegate source
  ParametersType delegateParameters = m_DelegateImageSource->GetParameters();
  unsigned int nextDelegateParameter = 0;
  for ( unsigned int i = 0; i < m_DelegateImageSource->GetNumberOfParameters(); ++i )
    {
    if ( m_EnabledArray[i] )
      {
      delegateParameters[i] = parameters[nextDelegateParameter++];
      }
    }

  if ( delegateParameters != m_DelegateImageSource->GetParameters() )
    {
    m_DelegateImageSource->SetParameters( delegateParameters );

    this->Modified();
    }
}

template< class TOutputImage >
typename MaskedParametricImageSource< TOutputImage >::ParametersType
MaskedParametricImageSource< TOutputImage >
::GetParameters() const
{
  if ( m_DelegateImageSource.GetPointer() == NULL )
    {
    itkExceptionMacro( << "Cannot get parameters when DelegateImageSource is not set" );
    }

  // Get only the active parameters in the delegate source
  ParametersType parameters( this->GetNumberOfParameters() );
  ParametersType delegateParameters = m_DelegateImageSource->GetParameters();
  unsigned int nextDelegateParameter = 0;
  for ( unsigned int i = 0; i < m_DelegateImageSource->GetNumberOfParameters(); ++i )
    {
    if ( m_EnabledArray[i] )
      {
      parameters[nextDelegateParameter++] = delegateParameters[i];
      }
    }

  return parameters;
}

template< class TOutputImage >
unsigned int
MaskedParametricImageSource< TOutputImage >
::GetNumberOfParameters() const
{
  unsigned int numberOfParameters = 0;
  for ( size_t i = 0; i < m_EnabledArray.size(); ++i )
    {
    if ( m_EnabledArray[i] )
      {
      numberOfParameters++;
      }
    }

  return numberOfParameters;
}

template< class TOutputImage >
unsigned int
MaskedParametricImageSource< TOutputImage >
::GetDelegateNumberOfParameters() const
{
  if ( m_DelegateImageSource.GetPointer() == NULL )
    {
    itkExceptionMacro( << "DelegateImageSource not set" );
    }

  return m_DelegateImageSource->GetNumberOfParameters();
}

template< class TOutputImage >
void
MaskedParametricImageSource< TOutputImage >
::SetParameterEnabled(unsigned int parameterIndex, bool enabled)
{
  if ( m_DelegateImageSource.GetPointer() == NULL )
    {
    itkExceptionMacro( << "DelegateImageSource not set" );
    }

  if ( parameterIndex >= m_DelegateImageSource->GetNumberOfParameters() )
    {
    itkExceptionMacro( << "Parameter index is out of bounds" );
    }

  m_EnabledArray[parameterIndex] = enabled;
}

template< class TOutputImage >
bool
MaskedParametricImageSource< TOutputImage >
::GetParameterEnabled(unsigned int parameterIndex) const
{
  if ( m_DelegateImageSource.GetPointer() == NULL )
    {
    itkExceptionMacro( << "DelegateImageSource not set" );
    }

  if ( parameterIndex >= m_DelegateImageSource->GetNumberOfParameters() )
    {
    itkExceptionMacro( << "Parameter index is out of bounds" );
    }

  return m_EnabledArray[parameterIndex];
}

template< class TOutputImage >
void
MaskedParametricImageSource< TOutputImage >
::GenerateData()
{
  if ( m_DelegateImageSource.GetPointer() == NULL )
    {
    itkExceptionMacro( << "DelegateImageSource not set" );
    }

  m_DelegateImageSource->SetSize( this->GetSize() );
  m_DelegateImageSource->SetSpacing( this->GetSpacing() );
  m_DelegateImageSource->SetOrigin( this->GetOrigin() );
  m_DelegateImageSource->SetDirection( this->GetDirection() );

  m_DelegateImageSource->GraftOutput( this->GetOutput() );
  m_DelegateImageSource->Update();
  this->GraftOutput( m_DelegateImageSource->GetOutput() );
}

template< class TOutputImage >
void
MaskedParametricImageSource< TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  this->Superclass::PrintSelf( os, indent );
  os << indent << "DelegateImageSource: " << m_DelegateImageSource.GetPointer() << std::endl;
  os << indent << "EnabledArray: [";
  for ( size_t i = 0; i < m_EnabledArray.size()-1; ++i )
    {
    os << m_EnabledArray[i] << ", ";
    }
  os << m_EnabledArray[m_EnabledArray.size()-1] << "]" << std::endl;
}

} // end namespace itk

#endif
